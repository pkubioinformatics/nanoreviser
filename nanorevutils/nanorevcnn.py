# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - nanorevcnn
 
 @Time: 2020/4/14 5:11 PM
 
 @Author: lotuswang
 
 
"""

import warnings
warnings.filterwarnings('ignore')
from keras.layers import  BatchNormalization, Add, TimeDistributed, Conv1D


def Conv1d_BN(x, nb_filter, kernel_size, strides=1, padding='same', name=None):
    if name is not None:
        bn_name = name + '_bn'
        conv_name = name + '_conv'
    else:
        bn_name = None
        conv_name = None
    x = TimeDistributed(Conv1D(nb_filter, kernel_size, padding=padding, strides=strides, activation='relu', name=conv_name))(x)
    x = TimeDistributed(BatchNormalization(name=bn_name))(x)
    return x


def identity_Block(inpt, nb_filter, kernel_size, strides=1, with_conv_shortcut=False):
    x = Conv1d_BN(inpt, nb_filter=nb_filter, kernel_size=kernel_size, padding='same')
    x = Conv1d_BN(x, nb_filter=nb_filter, kernel_size=kernel_size, padding='same')
    if with_conv_shortcut:
        shortcut = Conv1d_BN(inpt, nb_filter=nb_filter, strides=strides, kernel_size=kernel_size)
        x = Add()([x, shortcut])
        return x
    else:
        x = Add()([x, inpt])
        return x

