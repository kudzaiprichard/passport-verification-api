"""
Text DNN is the DNN based on predefined embedding word vectors.
"""
from gensim.models import KeyedVectors
from pathlib import Path
import pandas as pd
import spacy
from ..utils import nlp
import numpy as np
from keras import utils
from time import time


DEFAULT_LANGUAGE_MODEL = "en_core_web_sm"


class TextDNNDataGenerator(utils.Sequence):
    """
    Read data from given file, and generate training data
    File format:
    - csv format, separated by comma
    - 1st column is the text
    - 2nd column is class labels, if multiple label, using space as separator
    """
 
    def __init__(self, datafile:Path, embeddings, batch_size, n_classes, random_seed=None, shuffle=True, language_model=DEFAULT_LANGUAGE_MODEL):
    
        self.shuffle = shuffle
        self.batch_size = batch_size
        self.n_classes = n_classes
        self.random_seed = random_seed or int(time())

        np.random.seed(self.random_seed)
        self._prepare_embeddings(embeddings)
        self._prepare_data(datafile, language_model)

        self._prepare_nlps(datafile, language_model)
        self._prepare_trainset()
        self.on_epoch_end()

    def _prepare_embeddings(self, embeddings):
        if isinstance(embeddings, KeyedVectors):
            self.embeddings = embeddings
        elif isinstance(embeddings, Path):
            self.embeddings = KeyedVectors.load_word2vec_format(embeddings.as_posix())
        else:
            raise AssertionError(f'embeddings should be KeyedVectors or Path')

    def _prepare_data(self, datafile, language_model):
        df = pd.read_csv(datafile)

        classes = set(map(str.strip, ' '.join(df.iloc[:, 1].tolist()).split(' ')))
        assert len(classes) <= self.n_classes, f'given n_classess={self.n_classes} while NO. of classes={len(classes)}.'
        self.classes = list(classes)
        
        doc2vec = lambda x: np.mean([self.embeddings[t.lemma_] for t in nlp.tokenize(x)], axis=0)
        to_index = lambda x: [self.classes.index(c.strip()) for c in x.split(' ')]
        onehot = lambda x: np.bitwise_or.reduce(utils.to_categorical(to_index(x), num_classes=self.n_classes, dtype='int32'), axis=0)

        df['vectors'] = df.iloc[:, 0].map(doc2vec)
        df['labels'] = df.iloc[:, 1].map(onehot)

        self.data = df[['vectors', 'labels']]
    
    def __len__(self):
        return self.data.shape[0]

    def __getitem__(self, index):

        istart, iend = index * self.batch_size, (index + 1) * self.batch_size
        
        return self.data.iloc[istart:iend, :]

        '''
        X = np.empty((len(batch_instance_ids), self.n_vocab_size))
        Y = np.empty((len(batch_instance_ids), self.n_vocab_size), dtype=int)
        for i, instance_id in enumerate(batch_instance_ids):
            X[i, ] = self._oh_encode_tokens([self.instances[i][0]])
            Y[i, ] = self._oh_encode_tokens(self.instances[i][1])

        return X, Y
        '''

    def on_epoch_end(self):
        if self.shuffle:
            self.data = self.data.sample(1)
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
