from keras.layers import Input, Dense
from keras import Model
import logging


logger = logging.getLogger(__name__)


def build(input_shape, output_shape, dtype, activation=None, name='linear_model'):
    x = Input(shape=input_shape, dtype=dtype, name="linear_input")
    y = Dense(output_shape, activation=activation, name="linear_dense")(x)
    return Model(input=x, output=y, name=name)


