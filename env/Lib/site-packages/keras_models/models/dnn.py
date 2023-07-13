from . import LinearModel
from keras.layers import Input, Dense, Dropout
from keras.models import Model, Sequential
from itertools import chain
import logging


logger = logging.getLogger(__name__)


def build(input_shape, output_shape, dtype, output_activation="softmax", dense_layers=[1024, 512], r_drop=0.5, name="DNN"):
    assert len(dense_layers) >= 1, 'dense layers should not be empty!'

    x = Input(shape=input_shape, dtype=dtype, name="deep_input")

    h = x
    for i, v in enumerate(dense_layers):
        h = Dense(dense_layers[i], activation='relu', name=f"deep_dense_{i}")(h)
        h = Dropout(r_drop, name=f"deep_dropout_{i}")(h)

    y = Dense(output_shape, activation=output_activation, name="deep_output")(h)
    return Model(x, y, name=name)
    



