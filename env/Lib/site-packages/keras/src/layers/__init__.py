# Copyright 2016 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================
"""Keras layers API."""

# isort: off
import tensorflow.compat.v2 as tf

from keras.src.engine.base_layer import Layer
from keras.src.engine.base_preprocessing_layer import PreprocessingLayer

# Generic layers.
from keras.src.engine.input_layer import Input
from keras.src.engine.input_layer import InputLayer
from keras.src.engine.input_spec import InputSpec
from keras.src.layers.activation.elu import ELU
from keras.src.layers.activation.leaky_relu import LeakyReLU
from keras.src.layers.activation.prelu import PReLU

# Activations layers.
from keras.src.layers.activation.relu import ReLU
from keras.src.layers.activation.softmax import Softmax
from keras.src.layers.activation.thresholded_relu import ThresholdedReLU
from keras.src.layers.attention.additive_attention import AdditiveAttention
from keras.src.layers.attention.attention import Attention

# Attention layers.
from keras.src.layers.attention.multi_head_attention import MultiHeadAttention

# Convolution layer aliases.
# Convolution layers.
from keras.src.layers.convolutional.conv1d import Conv1D
from keras.src.layers.convolutional.conv1d import Convolution1D
from keras.src.layers.convolutional.conv1d_transpose import Conv1DTranspose
from keras.src.layers.convolutional.conv1d_transpose import Convolution1DTranspose
from keras.src.layers.convolutional.conv2d import Conv2D
from keras.src.layers.convolutional.conv2d import Convolution2D
from keras.src.layers.convolutional.conv2d_transpose import Conv2DTranspose
from keras.src.layers.convolutional.conv2d_transpose import Convolution2DTranspose
from keras.src.layers.convolutional.conv3d import Conv3D
from keras.src.layers.convolutional.conv3d import Convolution3D
from keras.src.layers.convolutional.conv3d_transpose import Conv3DTranspose
from keras.src.layers.convolutional.conv3d_transpose import Convolution3DTranspose
from keras.src.layers.convolutional.depthwise_conv1d import DepthwiseConv1D
from keras.src.layers.convolutional.depthwise_conv2d import DepthwiseConv2D
from keras.src.layers.convolutional.separable_conv1d import SeparableConv1D
from keras.src.layers.convolutional.separable_conv1d import SeparableConvolution1D
from keras.src.layers.convolutional.separable_conv2d import SeparableConv2D
from keras.src.layers.convolutional.separable_conv2d import SeparableConvolution2D

# Core layers.
from keras.src.layers.core.activation import Activation
from keras.src.layers.core.dense import Dense
from keras.src.layers.core.einsum_dense import EinsumDense
from keras.src.layers.core.embedding import Embedding
from keras.src.layers.core.identity import Identity
from keras.src.layers.core.lambda_layer import Lambda
from keras.src.layers.core.masking import Masking
from keras.src.layers.core.tf_op_layer import ClassMethod
from keras.src.layers.core.tf_op_layer import InstanceMethod
from keras.src.layers.core.tf_op_layer import InstanceProperty
from keras.src.layers.core.tf_op_layer import SlicingOpLambda
from keras.src.layers.core.tf_op_layer import TFOpLambda

# Locally-connected layers.
from keras.src.layers.locally_connected.locally_connected1d import (
    LocallyConnected1D,
)
from keras.src.layers.locally_connected.locally_connected2d import (
    LocallyConnected2D,
)

# Merging functions.
# Merging layers.
from keras.src.layers.merging.add import Add
from keras.src.layers.merging.add import add
from keras.src.layers.merging.average import Average
from keras.src.layers.merging.average import average
from keras.src.layers.merging.concatenate import Concatenate
from keras.src.layers.merging.concatenate import concatenate
from keras.src.layers.merging.dot import Dot
from keras.src.layers.merging.dot import dot
from keras.src.layers.merging.maximum import Maximum
from keras.src.layers.merging.maximum import maximum
from keras.src.layers.merging.minimum import Minimum
from keras.src.layers.merging.minimum import minimum
from keras.src.layers.merging.multiply import Multiply
from keras.src.layers.merging.multiply import multiply
from keras.src.layers.merging.subtract import Subtract
from keras.src.layers.merging.subtract import subtract
from keras.src.layers.normalization.batch_normalization import (
    SyncBatchNormalization,
)

# Normalization layers.
from keras.src.layers.normalization.group_normalization import GroupNormalization
from keras.src.layers.normalization.layer_normalization import LayerNormalization
from keras.src.layers.normalization.unit_normalization import UnitNormalization
from keras.src.layers.normalization.spectral_normalization import (
    SpectralNormalization,
)  # noqa: E501

# Preprocessing layers.
from keras.src.layers.preprocessing.category_encoding import CategoryEncoding
from keras.src.layers.preprocessing.discretization import Discretization
from keras.src.layers.preprocessing.hashed_crossing import HashedCrossing
from keras.src.layers.preprocessing.hashing import Hashing

# Image preprocessing layers.
from keras.src.layers.preprocessing.image_preprocessing import CenterCrop
from keras.src.layers.preprocessing.image_preprocessing import RandomBrightness
from keras.src.layers.preprocessing.image_preprocessing import RandomContrast
from keras.src.layers.preprocessing.image_preprocessing import RandomCrop
from keras.src.layers.preprocessing.image_preprocessing import RandomFlip
from keras.src.layers.preprocessing.image_preprocessing import RandomHeight
from keras.src.layers.preprocessing.image_preprocessing import RandomRotation
from keras.src.layers.preprocessing.image_preprocessing import RandomTranslation
from keras.src.layers.preprocessing.image_preprocessing import RandomWidth
from keras.src.layers.preprocessing.image_preprocessing import RandomZoom
from keras.src.layers.preprocessing.image_preprocessing import Rescaling
from keras.src.layers.preprocessing.image_preprocessing import Resizing
from keras.src.layers.preprocessing.integer_lookup import IntegerLookup
from keras.src.layers.preprocessing.normalization import Normalization
from keras.src.layers.preprocessing.string_lookup import StringLookup
from keras.src.layers.preprocessing.text_vectorization import TextVectorization
from keras.src.layers.regularization.activity_regularization import (
    ActivityRegularization,
)
from keras.src.layers.regularization.alpha_dropout import AlphaDropout

# Regularization layers.
from keras.src.layers.regularization.dropout import Dropout
from keras.src.layers.regularization.gaussian_dropout import GaussianDropout
from keras.src.layers.regularization.gaussian_noise import GaussianNoise
from keras.src.layers.regularization.spatial_dropout1d import SpatialDropout1D
from keras.src.layers.regularization.spatial_dropout2d import SpatialDropout2D
from keras.src.layers.regularization.spatial_dropout3d import SpatialDropout3D

# Reshaping layers.
from keras.src.layers.reshaping.cropping1d import Cropping1D
from keras.src.layers.reshaping.cropping2d import Cropping2D
from keras.src.layers.reshaping.cropping3d import Cropping3D
from keras.src.layers.reshaping.flatten import Flatten
from keras.src.layers.reshaping.permute import Permute
from keras.src.layers.reshaping.repeat_vector import RepeatVector
from keras.src.layers.reshaping.reshape import Reshape
from keras.src.layers.reshaping.up_sampling1d import UpSampling1D
from keras.src.layers.reshaping.up_sampling2d import UpSampling2D
from keras.src.layers.reshaping.up_sampling3d import UpSampling3D
from keras.src.layers.reshaping.zero_padding1d import ZeroPadding1D
from keras.src.layers.reshaping.zero_padding2d import ZeroPadding2D
from keras.src.layers.reshaping.zero_padding3d import ZeroPadding3D

if tf.__internal__.tf2.enabled():
    from keras.src.layers.normalization.batch_normalization import (
        BatchNormalization,
    )
    from keras.src.layers.normalization.batch_normalization_v1 import (
        BatchNormalization as BatchNormalizationV1,
    )

    BatchNormalizationV2 = BatchNormalization
else:
    from keras.src.layers.normalization.batch_normalization import (
        BatchNormalization as BatchNormalizationV2,
    )
    from keras.src.layers.normalization.batch_normalization_v1 import (
        BatchNormalization,
    )

    BatchNormalizationV1 = BatchNormalization

# Kernelized layers.
from keras.src.layers.kernelized import RandomFourierFeatures

# Pooling layer aliases.
# Pooling layers.
from keras.src.layers.pooling.average_pooling1d import AveragePooling1D
from keras.src.layers.pooling.average_pooling1d import AvgPool1D
from keras.src.layers.pooling.average_pooling2d import AveragePooling2D
from keras.src.layers.pooling.average_pooling2d import AvgPool2D
from keras.src.layers.pooling.average_pooling3d import AveragePooling3D
from keras.src.layers.pooling.average_pooling3d import AvgPool3D
from keras.src.layers.pooling.global_average_pooling1d import GlobalAveragePooling1D
from keras.src.layers.pooling.global_average_pooling1d import GlobalAvgPool1D
from keras.src.layers.pooling.global_average_pooling2d import GlobalAveragePooling2D
from keras.src.layers.pooling.global_average_pooling2d import GlobalAvgPool2D
from keras.src.layers.pooling.global_average_pooling3d import GlobalAveragePooling3D
from keras.src.layers.pooling.global_average_pooling3d import GlobalAvgPool3D
from keras.src.layers.pooling.global_max_pooling1d import GlobalMaxPool1D
from keras.src.layers.pooling.global_max_pooling1d import GlobalMaxPooling1D
from keras.src.layers.pooling.global_max_pooling2d import GlobalMaxPool2D
from keras.src.layers.pooling.global_max_pooling2d import GlobalMaxPooling2D
from keras.src.layers.pooling.global_max_pooling3d import GlobalMaxPool3D
from keras.src.layers.pooling.global_max_pooling3d import GlobalMaxPooling3D
from keras.src.layers.pooling.max_pooling1d import MaxPool1D
from keras.src.layers.pooling.max_pooling1d import MaxPooling1D
from keras.src.layers.pooling.max_pooling2d import MaxPool2D
from keras.src.layers.pooling.max_pooling2d import MaxPooling2D
from keras.src.layers.pooling.max_pooling3d import MaxPool3D
from keras.src.layers.pooling.max_pooling3d import MaxPooling3D
from keras.src.layers.rnn.abstract_rnn_cell import AbstractRNNCell

# Recurrent layers.
from keras.src.layers.rnn.base_rnn import RNN
from keras.src.layers.rnn.simple_rnn import SimpleRNN
from keras.src.layers.rnn.simple_rnn import SimpleRNNCell
from keras.src.layers.rnn.stacked_rnn_cells import StackedRNNCells

if tf.__internal__.tf2.enabled():
    from keras.src.layers.rnn.gru import GRU
    from keras.src.layers.rnn.gru import GRUCell
    from keras.src.layers.rnn.gru_v1 import GRU as GRUV1
    from keras.src.layers.rnn.gru_v1 import GRUCell as GRUCellV1
    from keras.src.layers.rnn.lstm import LSTM
    from keras.src.layers.rnn.lstm import LSTMCell
    from keras.src.layers.rnn.lstm_v1 import LSTM as LSTMV1
    from keras.src.layers.rnn.lstm_v1 import LSTMCell as LSTMCellV1

    GRUV2 = GRU
    GRUCellV2 = GRUCell
    LSTMV2 = LSTM
    LSTMCellV2 = LSTMCell
else:
    from keras.src.layers.rnn.gru import GRU as GRUV2
    from keras.src.layers.rnn.gru import GRUCell as GRUCellV2
    from keras.src.layers.rnn.gru_v1 import GRU
    from keras.src.layers.rnn.gru_v1 import GRUCell
    from keras.src.layers.rnn.lstm import LSTM as LSTMV2
    from keras.src.layers.rnn.lstm import LSTMCell as LSTMCellV2
    from keras.src.layers.rnn.lstm_v1 import LSTM
    from keras.src.layers.rnn.lstm_v1 import LSTMCell

    GRUV1 = GRU
    GRUCellV1 = GRUCell
    LSTMV1 = LSTM
    LSTMCellV1 = LSTMCell

# Serialization functions.
from keras.src.layers import serialization

# Wrapper functions.
from keras.src.layers.rnn.base_wrapper import Wrapper
from keras.src.layers.rnn.bidirectional import Bidirectional

# RNN Cell wrappers.
from keras.src.layers.rnn.cell_wrappers import DeviceWrapper
from keras.src.layers.rnn.cell_wrappers import DropoutWrapper
from keras.src.layers.rnn.cell_wrappers import ResidualWrapper

# Convolutional-recurrent layers.
from keras.src.layers.rnn.conv_lstm1d import ConvLSTM1D
from keras.src.layers.rnn.conv_lstm2d import ConvLSTM2D
from keras.src.layers.rnn.conv_lstm3d import ConvLSTM3D
from keras.src.layers.rnn.cudnn_gru import CuDNNGRU

# cuDNN recurrent layers.
from keras.src.layers.rnn.cudnn_lstm import CuDNNLSTM
from keras.src.layers.rnn.time_distributed import TimeDistributed
from keras.src.layers.serialization import deserialize
from keras.src.layers.serialization import deserialize_from_json
from keras.src.layers.serialization import get_builtin_layer
from keras.src.layers.serialization import serialize


class VersionAwareLayers:
    """Utility to be used internally to access layers in a V1/V2-aware fashion.

    When using layers within the Keras codebase, under the constraint that
    e.g. `layers.BatchNormalization` should be the `BatchNormalization` version
    corresponding to the current runtime (TF1 or TF2), do not simply access
    `layers.BatchNormalization` since it would ignore e.g. an early
    `compat.v2.disable_v2_behavior()` call. Instead, use an instance
    of `VersionAwareLayers` (which you can use just like the `layers` module).
    """

    def __getattr__(self, name):
        serialization.populate_deserializable_objects()
        if name in serialization.LOCAL.ALL_OBJECTS:
            return serialization.LOCAL.ALL_OBJECTS[name]
        return super().__getattr__(name)

