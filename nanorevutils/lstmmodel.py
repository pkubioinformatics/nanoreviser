# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - lstmmodel
 
 @Time: 2019/10/7 10:41 AM
 
 @Author: lotuswang
 
 
"""

import keras
from keras.preprocessing.sequence import pad_sequences
from keras.models import Sequential, Model
from keras.layers import Input, Embedding, Activation, BatchNormalization, Dropout
from keras.layers import Dense, Lambda, Flatten, RepeatVector,Permute,Multiply, Dot, concatenate,Add
from keras.layers import Bidirectional,LSTM
from keras.layers import TimeDistributed, Conv1D, MaxPool1D
from keras.optimizers import RMSprop, Adam
from keras.initializers import RandomNormal, Identity
from keras.utils import np_utils
from keras.utils import plot_model
import keras.backend as K
from nanorevutils.nanorevcnn import identity_Block


NB_CLASS = 6
SENT_LEN = 11
SIGNEL_LEN = 50
VEC_LEN = 6


def get_model1():
    signal_input = Input(shape=(SENT_LEN, 50, 1), dtype='float', name='signal_input')

    identity1 = identity_Block(signal_input, 8, 3)
    identity1 = Dropout(0.2)(identity1)
    #     identity2=identity_Block(identity1,8,3)
    #     identity2= Dropout(0.2)(identity2)

    signal_x_out = TimeDistributed(Flatten())(identity1)
    signal_x_out = TimeDistributed(Dense(64, name='signal_x_out'))(signal_x_out)

    read_input = Input(shape=(SENT_LEN, 6), dtype='float', name='read_input')
    read_rnn1 = Bidirectional(LSTM(16, return_sequences=True, name='read_rnn1', activation='tanh'))(read_input)
    read_rnn1 = BatchNormalization()(read_rnn1)
    read_rnn2 = Bidirectional(LSTM(64, return_sequences=True, name='read_rnn11', activation='tanh'))(read_rnn1)
    read_rnn2 = BatchNormalization()(read_rnn2)
    total_input = concatenate([read_rnn2, signal_x_out], axis=-1)
    total_rnn1 = Bidirectional(LSTM(128, return_sequences=True, name='total_rnn1', activation='tanh'))(total_input)
    total_rnn1 = BatchNormalization()(total_rnn1)
    total_rnn2 = Bidirectional(LSTM(64, return_sequences=True, name='total_rnn2', activation='tanh'))(total_rnn1)
    # total_rnn2 = BatchNormalization()(total_rnn2)

    # attention
    #     a = Dense(6, activation='softmax', name='activation')(total_rnn2)
    total_output = Dense(128, activation='relu')(total_rnn2)
    total_output = Dense(32, activation='relu')(total_output)
    # feature = Dense(16, activation='relu', name='feature')(total_output)
    main_out = Dense(6, activation='relu', name='main_out')(total_output)
    #     attention_out = Dot(axes=-1,normalize='L2', name='attention_out')([main_out, a])
    flat_out = Flatten()(main_out)
    feature = Dense(16, activation='relu', name='feature')(flat_out)
    predict = Dense(6, activation='softmax', name='final_out')(feature)

    input_target = Input(shape=(1,))
    centers = Embedding(6, 16)(input_target)
    l2_loss = Lambda(lambda x: K.sum(K.square(x[0] - x[1][:, 0]), 1, keepdims=True), name='l2_loss1')(
        [feature, centers])

    att_model_train = Model(inputs=[signal_input, read_input, input_target],
                            outputs=[predict, l2_loss])
    att_model_train.compile(optimizer='adam',
                            loss=['sparse_categorical_crossentropy', lambda y_true, y_pred: y_pred],
                            loss_weights=[1., 0.4], metrics=['accuracy', ])

    att_model_predict = Model(inputs=[signal_input, read_input, ],
                              outputs=[predict, ])
    att_model_predict.compile(optimizer='adam',
                              loss=['sparse_categorical_crossentropy'],
                              metrics=['accuracy'])
    return att_model_predict


def get_model2():
    signal_input = Input(shape=(SENT_LEN, 50, 1), dtype='float', name='signal_input')

    identity1 = identity_Block(signal_input, 8, 3)
    identity1 = Dropout(0.2)(identity1)
    #     identity2=identity_Block(identity1,8,3)
    #     identity2= Dropout(0.2)(identity2)

    signal_x_out = TimeDistributed(Flatten())(identity1)
    signal_x_out = TimeDistributed(Dense(64, name='signal_x_out'))(signal_x_out)

    read_input = Input(shape=(SENT_LEN, 6), dtype='float', name='read_input')
    read_rnn1 = Bidirectional(LSTM(16, return_sequences=True, name='read_rnn1', activation='tanh'))(read_input)
    read_rnn1 = BatchNormalization()(read_rnn1)
    read_rnn2 = Bidirectional(LSTM(64, return_sequences=True, name='read_rnn11', activation='tanh'))(read_rnn1)
    read_rnn2 = BatchNormalization()(read_rnn2)
    total_input = concatenate([read_rnn2, signal_x_out], axis=-1)
    total_rnn1 = Bidirectional(LSTM(128, return_sequences=True, name='total_rnn1', activation='tanh'))(total_input)
    total_rnn1 = BatchNormalization()(total_rnn1)
    total_rnn2 = Bidirectional(LSTM(64, return_sequences=True, name='total_rnn2', activation='tanh'))(total_rnn1)
    # total_rnn2 = BatchNormalization()(total_rnn2)

    # attention
    #     a = Dense(6, activation='softmax', name='activation')(total_rnn2)
    total_output = Dense(128, activation='relu')(total_rnn2)
    total_output = Dense(32, activation='relu')(total_output)
    # feature = Dense(16, activation='relu', name='feature')(total_output)
    main_out = Dense(6, activation='relu', name='main_out')(total_output)
    #     attention_out = Dot(axes=-1,normalize='L2', name='attention_out')([main_out, a])
    flat_out = Flatten()(main_out)
    feature = Dense(16, activation='relu', name='feature')(flat_out)
    predict = Dense(5, activation='softmax', name='final_out')(feature)

    input_target = Input(shape=(1,))
    centers = Embedding(5, 16)(input_target)
    l2_loss = Lambda(lambda x: K.sum(K.square(x[0] - x[1][:, 0]), 1, keepdims=True), name='l2_loss1')(
        [feature, centers])

    att_model_train = Model(inputs=[signal_input, read_input, input_target],
                            outputs=[predict, l2_loss])
    att_model_train.compile(optimizer='adam',
                            loss=['sparse_categorical_crossentropy', lambda y_true, y_pred: y_pred],
                            loss_weights=[1., 0.4], metrics=['accuracy', ])

    att_model_predict = Model(inputs=[signal_input, read_input, ],
                              outputs=[predict, ])
    att_model_predict.compile(optimizer='adam',
                              loss=['sparse_categorical_crossentropy'],
                              metrics=['accuracy'])
    return att_model_predict


