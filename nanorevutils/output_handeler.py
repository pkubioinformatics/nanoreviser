# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - output_handeler
 
 @Time: 2019/4/14 4:46 PM
 
 @Author: lotuswang
 
 
"""
import os
from nanorevutils.nanorev_fast5_handeler import handle_opts, opt_main
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

label_to_base = {5:'A', 4:'G', 3:'T', 2:'C',1:'-', 0:'D'}

def prep_read_fasta(fast5_fn, read_fasta_fn, bases):
    """
    :param fast5_fn: ./input/fast5/id_98490_ch139_read1203_strand.fast5
    :param read_fasta_fn: ./input/tmp/id_98490_ch139_read1203_strand.fasta
    :param bases: ['G', 'T', 'T', 'G', 'C', 'T', 'T', 'C', 'G', 'T', 'T']
    :return:
    return None
    """
    try:
        fast5_fn = fast5_fn.split('/')[-1]
        # reads_fasta = ''
        reads_fasta = ">" + fast5_fn.replace(' ', '|||') + '\n' + \
                      ''.join(bases)
        # print(reads_fasta)
        read_fp = open(read_fasta_fn, 'w')
        read_fp.write(reads_fasta)
        read_fp.close()
    except Exception as e:
        raise NotImplementedError('Error in writing .fasta file')
    return True


def prep_read_fastq(fast5_fn, read_fastq_fn, bases, qul):
    try:
        # reads_fasta = ''
        fast5_fn = fast5_fn.split('/')[-1]
        reads_fastq = "@" + fast5_fn.replace(' ', '|||') + '\n' + \
                      ''.join(bases)  \
                      +'+\n' + \
                      ''.join(qul)
        # print(reads_fasta)
        read_fp = open(read_fastq_fn, 'w')
        read_fp.write(reads_fastq)
        read_fp.close()
    except Exception as e:
        raise NotImplementedError('Error in writing .fastq file')
    return True


def get_dna_qul(temp_dir):
    dna_seq = ''
    dna_qul = ''
    temp_dir = os.path.join(temp_dir,'workspace')
    for file in os.listdir(temp_dir):
        if file.endswith('.fastq'):
            file = os.path.join(temp_dir,file)
            # print(file)
            with open(file, 'r') as out_fp:
                out_fp.seek(0)
                fastq_output = out_fp.readlines()
                # print(len(fastq_output))
                dna_seq=fastq_output[1][13:-13]
                dna_qul=fastq_output[3][13:-13]
        else:
            continue
    return dna_seq, dna_qul


def get_base_1(event_bases, y_pre, y_pre2):
    result = list()
    y_pre2=y_pre2-1
    result.append(label_to_base[y_pre[0]])
    for y_tmp, y_tmp2, base in zip(y_pre, y_pre2, event_bases):
        y_tmp = label_to_base.get(y_tmp,0)
        y_tmp2 = label_to_base.get(y_tmp2,0)
        if y_tmp==y_tmp2 and y_tmp in ['A','T','C','G']:
            result.append(y_tmp)
        elif y_tmp=='D'and y_tmp2 in ['A','T','C','G']:
            result.append(base)
            result.append(y_tmp2)
        elif y_tmp=='-' and y_tmp2=='-':
            # result.append(base)
            continue
        else:
            result.append(base)
    result = ''.join([tmp for tmp in result if tmp!='-'])
    return result


def get_base_2(event_bases, y_pre, y_pre2):
    result = list()
    result.append(y_pre[0])
    for y_tmp, y_tmp2, base in zip(y_pre, y_pre2, event_bases):
#         y_tmp = label_to_base.get(y_tmp,0)
#         y_tmp2 = label_to_base.get(y_tmp2,0)
        if y_tmp==y_tmp2 and y_tmp in ['A','T','C','G']:
            result.append(y_tmp)
        elif y_tmp=='D'and y_tmp2 in ['A','T','C','G']:
            result.append(base)
            result.append(y_tmp2)
        elif y_tmp=='-' and y_tmp2=='-':
            # result.append(base)
            continue
        else:
            result.append(base)
    result = ''.join([tmp for tmp in result if tmp!='-'])
    return result


def get_base_l(default_path, fast5_fn, temp_dir, event_bases=0, y_pred=0, y_pred2=0):
    result_DNA = ''
    result_qulity = ''
    input_dir = temp_dir
    opts = handle_opts(default_path, input_dir, temp_dir)
    exitFlag = opt_main(opts)
    if exitFlag==0:
        result_DNA, result_qulity = get_dna_qul(temp_dir)
    else:
        raise NotImplementedError('Error in revising file, like a broken .fast5 file.')
    assert (result_DNA != '')
    # assert (len(result_DNA) == len(result_qulity))
    return result_DNA, result_qulity


NB_CLASS = 6
SENT_LEN = 13
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



