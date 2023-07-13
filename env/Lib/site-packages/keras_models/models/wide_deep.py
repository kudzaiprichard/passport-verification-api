from . import LinearModel, DNN
from keras.layers import Input, Dense, Dropout, Concatenate
from keras.models import Model
import logging


logger = logging.getLogger(__name__)


def build(input_shape, output_shape, dtype, dense_layers=[1024, 512, 128], r_drop=0.5, name="wide_deep"):
    assert len(dense_layers) >= 2, 'invalid dense layers'

    x = Input(shape=input_shape, dtype=dtype, name="deep_wide_input")

    deep_model = DNN(input_shape, dense_layers[-1], dtype, output_activation=None, dense_layers=dense_layers[:-1], r_drop=r_drop, name="deep")
    deep_model.compile(loss='mae', optimizer='adam')
    yd = deep_model(x)

    xw = Concatenate()([x, yd])
    wide_model = LinearModel((input_shape[-1] + dense_layers[-1], ), output_shape, dtype, name="wide")
    wide_model.compile(loss='mae', optimizer='adam')
    y = wide_model(xw)
    
    return Model(x, y)

    


