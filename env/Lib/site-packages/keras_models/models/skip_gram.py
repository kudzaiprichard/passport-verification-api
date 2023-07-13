from keras.layers import Dense, Input
from keras.models import Model
from pathlib import Path
import spacy
from itertools import product, chain
from keras.utils import Sequence
from collections import Counter
import numpy as np
from time import time
from keras import backend as K


class SkipGramDataGenerator(Sequence):
 
    def __init__(self, datafile:Path, batch_size, r_sample=0.001, n_window=2, language_model="xx_ent_wiki_sm", random_seed=None, shuffle=True):
    
        self.shuffle = shuffle
        self.n_window = n_window
        self.r_sample = r_sample
        self.batch_size = batch_size

        np.random.seed(random_seed or int(time()))
        self._prepare_nlps(datafile, language_model)
        self._prepare_trainset()
        self.on_epoch_end()

    def _prepare_nlps(self, datafile, language_model):
        nlp = spacy.load(language_model)
        nlp.add_pipe(nlp.create_pipe('sentencizer'))

        with open(datafile, 'r') as f:
            doc = nlp(f.read())

        tokens_cnter = Counter([t.text for t in doc])
        self.kv_token_freq = {t: v / sum(tokens_cnter.values()) for t, v in tokens_cnter.items()}
        self.tokens = list(tokens_cnter.keys())
        self.n_vocab_size = len(self.tokens)

        P_subsample = lambda token: (np.sqrt(self.kv_token_freq[token] / self.r_sample) + 1) * self.r_sample / self.kv_token_freq[token]
        kv_token__P_subsample = {t: P_subsample(t) for t in self.tokens}

        subsample_sentence = lambda sent: [t for t in sent if kv_token__P_subsample[t.text] >= np.random.rand()]
        self.doc = [subsample_sentence(sent) for sent in doc.sents]

    def _prepare_trainset(self):
    
        offsets = list(range(-1 * self.n_window, 0)) + list(range(1, self.n_window + 1))
        pair_locs = lambda loc_ceil: [(i, offset) for i, offset in product(range(loc_ceil), offsets) if i + offset >=0 and i + offset < loc_ceil]

        instances = dict()
        for sentence in self.doc:
            for i, offset in pair_locs(len(sentence)):
                centra = sentence[i].text
                surround = sentence[i + offset].text
                instances[centra] = instances.get(centra, []) + [surround]
        self.instances = list(instances.items())

        self.kv_surround_token__instance_id = dict()
        for i, (_, surrounds) in enumerate(self.instances):
            for s in surrounds:
                self.kv_surround_token__instance_id[s] = self.kv_surround_token__instance_id.get(s, []) + [i]
    
    def __len__(self):
        return len(self.instances)

    def __getitem__(self, index):

        pos_instance_ids = self._pick_instances(self.batch_size // 2)
        neg_instance_ids = self._negative_sampling(pos_instance_ids)
        batch_instance_ids = pos_instance_ids + neg_instance_ids
        if len(batch_instance_ids) < self.batch_size:
            batch_instance_ids += self._pick_instances(self.batch_size - len(batch_instance_ids))

        X = np.empty((len(batch_instance_ids), self.n_vocab_size))
        Y = np.empty((len(batch_instance_ids), self.n_vocab_size), dtype=int)
        for i, instance_id in enumerate(batch_instance_ids):
            X[i, ] = self._oh_encode_tokens([self.instances[i][0]])
            Y[i, ] = self._oh_encode_tokens(self.instances[i][1])

        return X, Y
    
    def _oh_encode_tokens(self, tokens):
        encoding = np.zeros(self.n_vocab_size)
        indexs = [self.tokens.index(t) for t in tokens]
        encoding[indexs] = 1
        return encoding

    def _pick_instances(self, size):
        unfetched_instance_ids = np.squeeze(self.__unfetched_instance_flags.nonzero())
        if self.shuffle:
            np.random.shuffle(unfetched_instance_ids)

        picked = unfetched_instance_ids[:size]
        self.__unfetched_instance_flags[picked] = 0
        return picked.tolist()

    def _negative_sampling(self, positive_instance_ids):
    
        pos_surrounds = chain(*[self.instances[i][1] for i in positive_instance_ids]) 

        neg_candidate_instance_ids = []
        for token in set(self.kv_surround_token__instance_id.keys()) - set(pos_surrounds):
            for i in self.kv_surround_token__instance_id[token]:
                if self.__unfetched_instance_flags[i]:
                    neg_candidate_instance_ids.append(i)
            
        if self.shuffle:
            np.random.shuffle(neg_candidate_instance_ids)

        candidates_ids = neg_candidate_instance_ids[:self.batch_size // 2]
        self.__unfetched_instance_flags[candidates_ids] = 0
        return candidates_ids

    def on_epoch_end(self):
        self.__unfetched_instance_flags = np.ones(len(self.instances))
        return super().on_epoch_end()

    def summary(self):
        return f'''>>> Skip Gram Data Generator Summary:
        | Name             |                       Value                        |
        |------------------+----------------------------------------------------|
        | n_vocab_size     | {self.n_vocab_size:^50d} |
        | batch_size       | {self.batch_size:^50d} |
        | r_sample         | {self.r_sample:^50f} |
        | n_window         | {self.n_window:^50d} |
        | steps_each_epoch | {self.n_vocab_size // self.batch_size + 1:^50d} |
        \n'''


class SkipGram(Model):

    def __init__(self, n_vocab_size, n_embedding):
        super(SkipGram, self).__init__(name='SkipGram')
        
        self.n_vocab_size = n_vocab_size
        self.n_embedding = n_embedding
        self._prepare_network()

    def _prepare_network(self):
        self.hidden_layer = Dense(self.n_embedding, name="hidden_layer")
        self.output_layer = Dense(self.n_vocab_size, name="output_layer")
        
    def call(self, inputs):
        return self.output_layer(self.hidden_layer(inputs))

    def get_embeddings(self, inputs):
        return  K.function([self.hidden_layer.input],[self.hidden_layer.output])([inputs])[0]
