import spacy
from itertools import filterfalse


DEFAULT_LANGUAGE_MODEL = 'en_core_web_sm'


def tokenize(text, language_model=DEFAULT_LANGUAGE_MODEL, lemmatize=True, remove_stopwords=True):
    nlp = spacy.load(language_model)
    
    is_stopword = lambda x: remove_stopwords and x.norm_ in spacy.lang.en.stop_words.STOP_WORDS
    is_useless = lambda x: x.is_punct or x.is_digit or x.is_space or x.like_num or x.like_url or x.like_email
    is_removable = lambda x: is_stopword(x) or is_useless(x)

    token_value = lambda x: x.lemma_ if lemmatize else x.text
    yield from map(token_value, filterfalse(is_removable, nlp(text)))


def sentencizer(text, language_model=DEFAULT_LANGUAGE_MODEL):
    nlp = spacy.load(language_model)
    nlp.add_pipe(nlp.create_pipe('sentencizer'), before="parser")

    yield from map(lambda x: x.text.strip(), nlp(text).sents)
