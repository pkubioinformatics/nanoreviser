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
from nanorevutils.fileoptions import check_path, copy_file
from nanorevutils.nanorev_fast5_handeler import get_read_data, extract_fastq
from nanorevutils.preprocessing import signal_segmentation, get_base_color, get_base_label
from nanorevutils.input_handeler import parse_fasta
from nanorevutils.lstmmodel import get_model1, get_model2


def get_args():
    optParser = OptionParser(usage="%prog [-d] [-o]", version="%prog 1.0",
                             description="An Error-correction Tool for Nanopore Sequencing Based on a Deep Learning Algorithm")
    optParser.add_option('-d', '--fast5_base_dir', action='store', type="string",
                         dest='fast5_base_dir',
                         help='path to the fast5 files')
    optParser.add_option('-o', '--output_dir', action='store', type="string", dest='output_dir',
                         default='./unitest/nanorev_training_result/',
                         help='path to store the output sumarry files')
    optParser.add_option('-r', '--reference', action='store', type="string", dest='genome_fn',
                         default='./unitest/nanorev_training_result/',
                         help='path to store the output files')
    optParser.add_option('-m', '--output_model', action='store', type="string", dest='model_dir',
                         help='path to store the output files')

    optParser.add_option('-L', '--output_format', action='store', type="string", dest='output_format',
                         default='sam',)
    optParser.add_option("--thread", action="store", type="int", dest="thread",
                         default=1,
                         help='thread, default is 1')
    optParser.add_option('-t', '--tmp_dir', action='store', type="string", dest='temp_dir',
                          default='./tmp/',
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


    optParser.add_option("-b", "--batch_size", action="store", type="int", dest="batch_size",
                         default=256,
                         help='batch size, default is 256')
    optParser.add_option("-e", "--epochs", action="store", type="int", dest="epochs",
                         default=50,
                         help='epochs, default is 50')
    optParser.add_option("-w", "--window_size", action="store", type="int", dest="window_size",
                         default=13,
                         help='window size, default is 13')
    optParser.add_option("-c", "--read_count", action="store", type="int", dest="read_count",
                         default=0,
                         help='the number of read included in the training data, must '
                              'smaller than the number of files stored in fast5_base_dir, '
                              '0 for use all the files in the fast5_base_dir and defult is 0.')
    optParser.add_option("--validation_split", action="store", type="float", dest="validation_split",
                         default=0.01,
                         help='validation data size, default is 0.01, which means 1% '
                              'reads in fast5_base_dir would be used as validation data.')


    optParser.add_option('--test_mode', action='store_true', default=False,
                        help='just for unitest')


    optParser.add_option('--model1_predict_dir', action='store', type="string", dest='model1_predict_dir',
                         default='./model/ecoli_win13_50ep_model1.h5',
                        help='model dirs for model1')
    optParser.add_option('--model2_predict_dir', action='store', type="string", dest='model2_predict_dir',
                         default='./model/ecoli_win13_50ep_model2.h5',
                        help='model dirs for model2')
    optParser.add_option("-v", "--virsion", action="store_true", dest="virsion",
                         help="version of NanoReviser")

    (tmp_args, _) = optParser.parse_args()
    if tmp_args.virsion:
        print("The virsion of NanoReviser : 0.1 ")
    return tmp_args


def this_folder():
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    else:
        return os.path.dirname(__file__)
default_path = get_default_path(this_folder())



if __name__ == '__main__':
    ar_args=get_args()
