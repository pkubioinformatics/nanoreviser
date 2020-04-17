# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - NanoReviser_train
 
 @Time: 2020/4/15 5:32 PM
 
 @Author: lotuswang
 
 
"""
import warnings

warnings.filterwarnings('ignore')

import os
import time
import sys
from optparse import OptionParser
from multiprocessing import Pool
import shutil
import numpy as np
import pandas as pd

from nanorevutils.nanolog import logger_config
from albacore.path_utils import get_default_path
from nanorevutils.fileoptions import check_path, summary_generate, write_sumery_file, model_fn_generate
from nanorevutils.lstmmodel import get_model1, get_model2
from nanorevutils.nanorevtrainutils import train_preprocessing, get_trainning_input

def get_args():
    optParser = OptionParser(usage="%prog [-d] [-o]", version="%prog 1.0",
                             description="An Training Tool for NanoReviser")
    optParser.add_option('-d', '--fast5_base_dir', action='store', type="string",
                         default='./unitest/training_data/fast5/',
                         dest='fast5_base_dir',
                         help='path to the fast5 files')
    optParser.add_option('-o', '--output_dir', action='store', type="string", dest='output_dir',
                         default='./unitest/nanorev_training_result/',
                         help='path to store the output sumarry files')
    optParser.add_option('-r', '--reference', action='store', type="string", dest='genome_fn',
                         default='./unitest/training_data/reference.fasta',
                         help='path to store the output files')

    optParser.add_option('--model_type', action='store', type="string", dest='model_type',
                         default='both',
                         help='both, model1 or model2, default is both')
    optParser.add_option('-S', '--species', action='store', type="string", dest='species',
                         default='unitest',
                         help='species of training data, default is unitest for test Nanoreviser_train')
    optParser.add_option('-M', '--output_model', action='store', type="string", dest='model_dir',
                         default='./model/',
                         help='path to store the model files')

    optParser.add_option('-m', '--mapper_exe', action='store', type="string", dest='graphmap_exe',
                         default='graphmap',
                         help='the align tool for generate the lable of training data, default is graphmap')
    optParser.add_option('-L', '--output_format', action='store', type="string", dest='output_format',
                         default='sam',)
    optParser.add_option("--thread", action="store", type="int", dest="thread",
                         default=1,
                         help='thread, default is 1')

    optParser.add_option('-t', '--tmp_dir', action='store', type="string", dest='temp_dir',
                          default='./train_tmp/',
                         help='path to the tmp dir, which is used to store the preprocessing files')
    optParser.add_option('-f', '--failed_read', action='store', type="string", dest='failed_reads_filename',
                          default='failed_reads.txt',
                         help='document to log the failed reads, default is failed_read.txt')
    optParser.add_option('-g', '--basecall_group', action='store', type="string", dest='basecall_group',
                         default='Basecall_1D_000',
                         help='attrs for finding the events file in fast5 file, default is Basecall_1D_000')
    optParser.add_option('-s', '--basecall_subgroup', action='store', type="string", dest='basecall_subgroup',
                         default='BaseCalled_template',
                         help='attrs for finding the events file in fast5 file, default is BaseCalled_template')

    #options for deeplearning
    optParser.add_option("-b", "--batch_size", action="store", type="int", dest="batch_size",
                         default=512,
                         help='batch size, default is 256')
    optParser.add_option("-e", "--epochs", action="store", type="int", dest="epochs",
                         default=2,
                         help='epochs, default is 50')
    optParser.add_option("-w", "--window_size", action="store", type="int", dest="window_size",
                         default=13,
                         help='window size, default is 13')
    optParser.add_option("-c", "--read_counts", action="store", type="int", dest="read_counts",
                         default=1,
                         help='the number of read included in the training data, must '
                              'smaller than the number of files stored in fast5_base_dir, '
                              '0 for use all the files in the fast5_base_dir and defult is 0.')
    optParser.add_option("--validation_split", action="store", type="float", dest="validation_split",
                         default=0.01,
                         help='validation data size, default is 0.01, which means 1% '
                              'reads in fast5_base_dir would be used as validation data.')

    optParser.add_option('--model1_train_dir', action='store', type="string", dest='model1_train_dir',
                         default='',
                        help='model dirs for trained model1, for transfer learning')
    optParser.add_option('--model2_train_dir', action='store', type="string", dest='model2_train_dir',
                         default='',
                        help='model dirs for trained model1, for transfer learning')

    optParser.add_option('--test_mode', action='store_true', default=False,
                         help='just for unitest')
    optParser.add_option("-v", "--virsion", action="store_true", dest="virsion",
                         help="version of NanoReviser")

    (tmp_args, _) = optParser.parse_args()
    tmp_args.model_dir = str(tmp_args.model_dir) +'/'+ str(tmp_args.species)+'/'
    tmp_args.train_input_dir = str(tmp_args.model_dir) + '/training_input/'
    tmp_args.train_model_dir = str(tmp_args.model_dir) + '/training_model/'
    if tmp_args.virsion:
        print("The virsion of NanoReviser : 1.0 ")
    return tmp_args


def this_folder():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    else:
        return os.path.dirname(__file__)
default_path = get_default_path(this_folder())


if __name__ == '__main__':
    ar_args=get_args()
    if ar_args.test_mode:
        logger = logger_config(log_path='./unitest/unitest_log.txt', logging_name='unitest')
        ar_args.epochs=2
        ar_args.read_counts=1
        ar_args.window_size=5
    try:
        start_time = time.time()
        try:
            shutil.rmtree(ar_args.temp_dir)
        except Exception:
            tttt = 1
        check_path(ar_args.temp_dir)
        check_path(ar_args.output_dir)
        check_path(ar_args.train_input_dir)
        train_preprocessing(ar_args)
        check_path(ar_args.train_model_dir)
        try:
            x_train, signal_x_train, y_train, y_train2 = get_trainning_input(ar_args.test_mode,
                                                                             ar_args.train_input_dir,
                                                                             ar_args.window_size)
        except Exception as e:
            raise RuntimeError(e)


        try:
            att_model_train1, att_model_predict1 = get_model1(ar_args.window_size)
            if not ar_args.test_mode:
                print('[p:::] model buiding.....')
        except Exception as e:
            raise RuntimeError('！！！[Error] loading training model ',e)
        try:
            if not ar_args.test_mode:
                print('[p:::] start to training model1, please waiting......')
                verbose=2
            else:
                verbose=0
            m1_start_t = time.time()
            att_model_train1, att_model_predict1 = get_model1(ar_args.window_size)
            history1 = att_model_train1.fit([signal_x_train[:, :, :, np.newaxis], x_train, y_train],
                                             [y_train, np.zeros((len(y_train2), 1))],
                                             class_weight={0: 3, 1: 5, 2: 1, 3: 1, 4: 1, 5: 1},
                                             validation_split=ar_args.validation_split,
                                             shuffle=True,
                                             epochs=ar_args.epochs,
                                             batch_size=ar_args.batch_size,
                                             verbose=verbose)

            model1_pre_fn, model1_train_fn, model1_history_fn, model1_summary_fn = model_fn_generate(ar_args,'model1')
            att_model_train1.save_weights(model1_train_fn)
            att_model_predict1.save_weights(model1_pre_fn)
            model1_summary = summary_generate(ar_args,m1_start_t)
            write_sumery_file(dict(history1.history), model1_summary, model1_history_fn, model1_summary_fn)
            if not ar_args.test_mode:
                print('[p:::] model 1 completed......')
        except Exception as e:
            raise RuntimeError('！！！[Error]training model 1...... ', e)

        try:
            if not ar_args.test_mode:
                print('[p:::] start to training model2, please waiting......')
                verbose = 2
            else:
                verbose = 0
            m2_start_t=time.time()
            att_model_train2, att_model_predict2 = get_model2(ar_args.window_size)
            history2 = att_model_train2.fit([signal_x_train[:, :, :, np.newaxis], x_train, y_train2],
                                            [y_train2, np.zeros((len(y_train2), 1))],
                                            class_weight={0: 3, 1: 5, 2: 1, 3: 1, 4: 1, 5: 1},
                                            validation_split=ar_args.validation_split,
                                            shuffle=True,
                                            epochs=ar_args.epochs,
                                            batch_size=ar_args.batch_size,
                                            verbose=verbose)

            model2_pre_fn, model2_train_fn, model2_history_fn, model2_summary_fn = model_fn_generate(ar_args,
                                                                                                     'model2')
            att_model_train2.save_weights(model2_train_fn)
            att_model_predict2.save_weights(model2_pre_fn)
            model2_summary = summary_generate(ar_args, m2_start_t)
            write_sumery_file(dict(history2.history), model2_summary, model2_history_fn, model2_summary_fn)
            if not ar_args.test_mode:
                print('[p:::] model 2 completed......')

        except Exception as e:
            raise RuntimeError('！！！[Error]training model 2', e)

        try:
            end_time = time.time()
            if not ar_args.test_mode:
                print('[s:::] The training time of NanoReviser_train is :%.2f seconds' % (end_time - start_time))
            else:
                logger.info("Congratulations, NanoReviser_train is installed properly")
                shutil.rmtree(ar_args.output_dir)
                shutil.rmtree(ar_args.model_dir)
            shutil.rmtree(ar_args.temp_dir)
        except Exception as e:
            print('！！！[Error] remove tmp dir ' + ar_args.temp_dir + e)
    except Exception as e:
        if ar_args.test_mode:
            logger.error(e)
        else:
            print(e)
