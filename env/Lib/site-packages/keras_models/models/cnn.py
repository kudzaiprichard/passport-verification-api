from keras.layers import Dense, Conv2D, MaxPool2D, Dropout, Flatten
import logging
from keras.models import Sequential


logger = logging.getLogger(__name__)


def build(input_shape, num_classes, filters, kernel_size, pool_size, padding, r_dropout):
    return Sequential([
        Conv2D(filters[0], kernel_size, padding=padding, input_shape=input_shape, activation='relu', name='conv1_1'),
        Conv2D(filters[0], kernel_size, activation='relu', name='conv1_2'),
        MaxPool2D(pool_size, name='pool1'),
        Dropout(r_dropout, name='dropout1'),
        Conv2D(filters[1], kernel_size, padding=padding, activation='relu', name='conv2_1'),
        Conv2D(filters[1], kernel_size, activation='relu', name='conv2_2'),
        MaxPool2D(pool_size, name='pool2'),
        Dropout(r_dropout, name='dropout2'),
        Flatten(),
        Dense(512, activation='relu'),
        Dropout(r_dropout, name='dropout3'),
        Dense(num_classes, activation='softmax')
    ])
