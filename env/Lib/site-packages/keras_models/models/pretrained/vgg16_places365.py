# -*- coding: utf-8 -*-
'''VGG16-places365 model for Keras
This code is forked from: https://github.com/GKalliatakis/Keras-VGG16-places365/blob/master/vgg16_places_365.py
To be re-wrote in Keras way in future.

# Reference:
- [Places: A 10 million Image Database for Scene Recognition](http://places2.csail.mit.edu/PAMI_places.pdf)
'''

import os

import warnings
import numpy as np
from PIL import Image
from cv2 import resize
from pathlib import Path

from keras import backend as K
from keras.layers import Input
from keras.layers.core import Activation, Dense, Flatten
from keras.layers.pooling import MaxPooling2D
from keras.models import Model
from keras.layers import Conv2D
from keras.regularizers import l2
from keras.layers.core import Dropout
from keras.layers import GlobalAveragePooling2D
from keras.layers import GlobalMaxPooling2D
from keras_applications.imagenet_utils import _obtain_input_shape
from keras.engine.topology import get_source_inputs
from keras.utils.data_utils import get_file
from keras.utils import layer_utils
from keras.preprocessing import image
from keras.applications.imagenet_utils import preprocess_input

WEIGHTS_PATH = 'https://github.com/Marcnuth/Keras-Models/releases/download/v0.0.5/vgg16-places365_weights_tf_dim_ordering_tf_kernels.h5'
CLASS_LABELS = [
    'airfield', 'airplane_cabin', 'airport_terminal', 'alcove', 'alley', 'amphitheater', 'amusement_arcade', 'amusement_park', 'apartment_building/outdoor', 
    'aquarium', 'aqueduct', 'arcade', 'arch', 'archaelogical_excavation', 'archive', 'arena/hockey', 'arena/performance', 'arena/rodeo', 'army_base', 
    'art_gallery', 'art_school', 'art_studio', 'artists_loft', 'assembly_line', 'athletic_field/outdoor', 'atrium/public', 'attic', 'auditorium', 
    'auto_factory', 'auto_showroom', 'badlands', 'bakery/shop', 'balcony/exterior', 'balcony/interior', 'ball_pit', 'ballroom', 'bamboo_forest', 'bank_vault', 
    'banquet_hall', 'bar', 'barn', 'barndoor', 'baseball_field', 'basement', 'basketball_court/indoor', 'bathroom', 'bazaar/indoor', 'bazaar/outdoor', 'beach', 
    'beach_house', 'beauty_salon', 'bedchamber', 'bedroom', 'beer_garden', 'beer_hall', 'berth', 'biology_laboratory', 'boardwalk', 'boat_deck', 'boathouse', 
    'bookstore', 'booth/indoor', 'botanical_garden', 'bow_window/indoor', 'bowling_alley', 'boxing_ring', 'bridge', 'building_facade', 'bullring', 
    'burial_chamber', 'bus_interior', 'bus_station/indoor', 'butchers_shop', 'butte', 'cabin/outdoor', 'cafeteria', 'campsite', 'campus', 'canal/natural', 
    'canal/urban', 'candy_store', 'canyon', 'car_interior', 'carrousel', 'castle', 'catacomb', 'cemetery', 'chalet', 'chemistry_lab', 'childs_room', 
    'church/indoor', 'church/outdoor', 'classroom', 'clean_room', 'cliff', 'closet', 'clothing_store', 'coast', 'cockpit', 'coffee_shop', 'computer_room', 
    'conference_center', 'conference_room', 'construction_site', 'corn_field', 'corral', 'corridor', 'cottage', 'courthouse', 'courtyard', 'creek', 'crevasse', 
    'crosswalk', 'dam', 'delicatessen', 'department_store', 'desert/sand', 'desert/vegetation', 'desert_road', 'diner/outdoor', 'dining_hall', 'dining_room', 
    'discotheque', 'doorway/outdoor', 'dorm_room', 'downtown', 'dressing_room', 'driveway', 'drugstore', 'elevator/door', 'elevator_lobby', 'elevator_shaft', 
    'embassy', 'engine_room', 'entrance_hall', 'escalator/indoor', 'excavation', 'fabric_store', 'farm', 'fastfood_restaurant', 'field/cultivated', 'field/wild', 
    'field_road', 'fire_escape', 'fire_station', 'fishpond', 'flea_market/indoor', 'florist_shop/indoor', 'food_court', 'football_field', 'forest/broadleaf', 
    'forest_path', 'forest_road', 'formal_garden', 'fountain', 'galley', 'garage/indoor', 'garage/outdoor', 'gas_station', 'gazebo/exterior', 
    'general_store/indoor', 'general_store/outdoor', 'gift_shop', 'glacier', 'golf_course', 'greenhouse/indoor', 'greenhouse/outdoor', 'grotto', 
    'gymnasium/indoor', 'hangar/indoor', 'hangar/outdoor', 'harbor', 'hardware_store', 'hayfield', 'heliport', 'highway', 'home_office', 'home_theater', 
    'hospital', 'hospital_room', 'hot_spring', 'hotel/outdoor', 'hotel_room', 'house', 'hunting_lodge/outdoor', 'ice_cream_parlor', 'ice_floe', 'ice_shelf', 
    'ice_skating_rink/indoor', 'ice_skating_rink/outdoor', 'iceberg', 'igloo', 'industrial_area', 'inn/outdoor', 'islet', 'jacuzzi/indoor', 'jail_cell', 
    'japanese_garden', 'jewelry_shop', 'junkyard', 'kasbah', 'kennel/outdoor', 'kindergarden_classroom', 'kitchen', 'lagoon', 'lake/natural', 'landfill', 
    'landing_deck', 'laundromat', 'lawn', 'lecture_room', 'legislative_chamber', 'library/indoor', 'library/outdoor', 'lighthouse', 'living_room', 'loading_dock', 
    'lobby', 'lock_chamber', 'locker_room', 'mansion', 'manufactured_home', 'market/indoor', 'market/outdoor', 'marsh', 'martial_arts_gym', 'mausoleum', 'medina', 
    'mezzanine', 'moat/water', 'mosque/outdoor', 'motel', 'mountain', 'mountain_path', 'mountain_snowy', 'movie_theater/indoor', 'museum/indoor', 'museum/outdoor', 
    'music_studio', 'natural_history_museum', 'nursery', 'nursing_home', 'oast_house', 'ocean', 'office', 'office_building', 'office_cubicles', 'oilrig', 
    'operating_room', 'orchard', 'orchestra_pit', 'pagoda', 'palace', 'pantry', 'park', 'parking_garage/indoor', 'parking_garage/outdoor', 'parking_lot', 
    'pasture', 'patio', 'pavilion', 'pet_shop', 'pharmacy', 'phone_booth', 'physics_laboratory', 'picnic_area', 'pier', 'pizzeria', 'playground', 'playroom', 
    'plaza', 'pond', 'porch', 'promenade', 'pub/indoor', 'racecourse', 'raceway', 'raft', 'railroad_track', 'rainforest', 'reception', 'recreation_room', 
    'repair_shop', 'residential_neighborhood', 'restaurant', 'restaurant_kitchen', 'restaurant_patio', 'rice_paddy', 'river', 'rock_arch', 'roof_garden', 
    'rope_bridge', 'ruin', 'runway', 'sandbox', 'sauna', 'schoolhouse', 'science_museum', 'server_room', 'shed', 'shoe_shop', 'shopfront', 'shopping_mall/indoor', 
    'shower', 'ski_resort', 'ski_slope', 'sky', 'skyscraper', 'slum', 'snowfield', 'soccer_field', 'stable', 'stadium/baseball', 'stadium/football', 
    'stadium/soccer', 'stage/indoor', 'stage/outdoor', 'staircase', 'storage_room', 'street', 'subway_station/platform', 'supermarket', 'sushi_bar', 'swamp', 
    'swimming_hole', 'swimming_pool/indoor', 'swimming_pool/outdoor', 'synagogue/outdoor', 'television_room', 'television_studio', 'temple/asia', 'throne_room', 
    'ticket_booth', 'topiary_garden', 'tower', 'toyshop', 'train_interior', 'train_station/platform', 'tree_farm', 'tree_house', 'trench', 'tundra', 
    'underwater/ocean_deep', 'utility_room', 'valley', 'vegetable_garden', 'veterinarians_office', 'viaduct', 'village', 'vineyard', 'volcano', 
    'volleyball_court/outdoor', 'waiting_room', 'water_park', 'water_tower', 'waterfall', 'watering_hole', 'wave', 'wet_bar', 'wheat_field', 'wind_farm', 
    'windmill', 'yard', 'youth_hostel', 'zen_garden']


def VGG16_Places365(model_file=None):
    """Instantiates the VGG16-places365 architecture.

    Optionally loads weights pre-trained on Places. Note that when using TensorFlow,
    for best performance you should set `image_data_format="channels_last"` in your Keras config at ~/.keras/keras.json.

    The model and the weights are compatible with both TensorFlow and Theano. 
    The data format convention used by the model is the one specified in your Keras config file.

    # Returns
        A Keras model instance.
    """

    # Determine proper input shape
    input_shape = _obtain_input_shape(None, default_size=224,  min_size=48, data_format=K.image_data_format(), require_flatten=True)
    img_input = Input(shape=input_shape)


    # Block 1
    x = Conv2D(filters=64, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block1_conv1')(img_input)
    x = Conv2D(filters=64, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block1_conv2')(x)
    x = MaxPooling2D(pool_size=(2, 2), strides=(2, 2), name="block1_pool", padding='valid')(x)

    # Block 2
    x = Conv2D(filters=128, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block2_conv1')(x)
    x = Conv2D(filters=128, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block2_conv2')(x)
    x = MaxPooling2D(pool_size=(2, 2), strides=(2, 2), name="block2_pool", padding='valid')(x)

    # Block 3
    x = Conv2D(filters=256, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block3_conv1')(x)
    x = Conv2D(filters=256, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block3_conv2')(x)
    x = Conv2D(filters=256, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block3_conv3')(x)
    x = MaxPooling2D(pool_size=(2, 2), strides=(2, 2), name="block3_pool", padding='valid')(x)

    # Block 4
    x = Conv2D(filters=512, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block4_conv1')(x)
    x = Conv2D(filters=512, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block4_conv2')(x)
    x = Conv2D(filters=512, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block4_conv3')(x)
    x = MaxPooling2D(pool_size=(2, 2), strides=(2, 2), name="block4_pool", padding='valid')(x)

    # Block 5
    x = Conv2D(filters=512, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block5_conv1')(x)
    x = Conv2D(filters=512, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block5_conv2')(x)
    x = Conv2D(filters=512, kernel_size=3, strides=(1, 1), padding='same', kernel_regularizer=l2(0.0002), activation='relu', name='block5_conv3')(x)
    x = MaxPooling2D(pool_size=(2, 2), strides=(2, 2), name="block5_pool", padding='valid')(x)

    # Classification block
    x = Flatten(name='flatten')(x)
    x = Dense(4096, activation='relu', name='fc1')(x)
    x = Dropout(0.5, name='drop_fc1')(x)

    x = Dense(4096, activation='relu', name='fc2')(x)
    x = Dropout(0.5, name='drop_fc2')(x)
    x = Dense(365, activation='softmax', name="predictions")(x)

    # Create model.
    model = Model(img_input, x, name='vgg16-places365')

    # load weights
    if model_file:
        model.load_weights(Path(model_file).absolute().as_posix())
    else:
        weights_path = get_file('vgg16-places365_weights_tf_dim_ordering_tf_kernels.h5', WEIGHTS_PATH, cache_subdir='models')
        model.load_weights(weights_path)

    if K.backend() == 'theano':
        layer_utils.convert_all_kernels_in_model(model)

    if K.image_data_format() == 'channels_first':
        maxpool = model.get_layer(name='block5_pool')
        shape = maxpool.output_shape[1:]
        dense = model.get_layer(name='fc1')
        layer_utils.convert_dense_weights_data_format(dense, shape, 'channels_first')

        if K.backend() == 'tensorflow':
            warnings.warn(
                'You are using the TensorFlow backend, yet you are using the Theano image data format convention (`image_data_format="channels_first"`). '
                'For best performance, set `image_data_format="channels_last"` in your Keras config at ~/.keras/keras.json.'
            )

    return model


def predict(image_files, n_top=5, model=None):

    model = model or VGG16_Places365()
    assert isinstance(image_files, list), 'image_files should be a list'
    assert len(image_files) > 0, 'image_files should not be empty'

    file2img = lambda f: np.expand_dims(resize(np.array(Image.open(f), dtype=np.uint8), (224, 224)), 0)
    images = [file2img(Path(f).absolute().as_posix()) for f in image_files]

    predict_scores = lambda img: np.argsort(model.predict(img)[0])[::-1][0:n_top]
    predict_labels = lambda img: [CLASS_LABELS[v] for v in predict_scores(img)]

    return [predict_labels(img) for img in images]



