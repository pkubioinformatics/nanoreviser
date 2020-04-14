# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - NanoReviser
 
 @Time: 2020/3/14 5:23 PM
 
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

from albacore.path_utils import get_default_path
from nanorevutils.fileoptions import check_path, copy_file
from nanorevutils.nanorev_fast5_handeler import get_read_data
from nanorevutils.preprocessing import signal_segmentation, get_base_color, get_base_label
from nanorevutils.input_handeler import parse_fasta
from nanorevutils.lstmmodel import get_model1, get_model2
from nanorevutils.output_handeler import get_base_l, prep_read_fasta, prep_read_fastq, get_base_1


os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


def get_args():
    optParser = OptionParser(usage="%prog [-d] [-o]", version="%prog 1.0")

    optParser.add_option('-d', '--fast5_base_dir', action='store', type="string",
                         dest='fast5_base_dir', )
    optParser.add_option('-o', '--output_dir', action='store', type="string", dest='output_dir',
                         default='./unitest/nanorev_output/')
    optParser.add_option('-F', '--output_format', action='store', type="string", dest='output_format',
                         default='fasta')
    optParser.add_option("--thread", action="store", type="int", dest="thread",
                         default=5)
    optParser.add_option('-m', '--mapper_dir', action='store', type="string", dest='graphmap_exe',
                         default='graphmap')
    optParser.add_option('-t', '--tmp_dir', action='store', type="string", dest='temp_dir',
                         default='./unitest/tmp/')
    optParser.add_option('-e', '--failed_read', action='store', type="string", dest='failed_reads_filename',
                         default='failed_reads.txt')
    optParser.add_option('-g', '--basecall_group', action='store', type="string", dest='basecall_group',
                         default='Basecall_1D_000')
    optParser.add_option('-s', '--basecall_subgroup', action='store', type="string", dest='basecall_subgroup',
                         default='BaseCalled_template')
    optParser.add_option('--model1_predict_dir', action='store', type="string", dest='model1_predict_dir',
                         default='./model/ecoli_win13_50ep_model1.h5')
    optParser.add_option('--model2_predict_dir', action='store', type="string", dest='model2_predict_dir',
                         default='./model/ecoli_win13_50ep_model2.h5')
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

def provide_fasta(name, fast5_fn_sg, args):

    print('[b:::] Run task %s pid %s' % (name, os.getpid()),' : ',fast5_fn_sg)
    fast5_fn = os.path.join(args.fast5_base_dir, fast5_fn_sg)
    basecall_tmp_dir = str(args.temp_dir) + str(name) + '/basecall_tmp/'
    copy_file(fast5_fn, basecall_tmp_dir)
    fast5_fn = os.path.join(basecall_tmp_dir, fast5_fn_sg)
    try:
        (read_start_rel_to_raw, starts_rel_to_read, event_length, event_bases,
         or_raw_signal, raw_mean, raw_std) = get_read_data(fast5_fn, args.basecall_group, args.basecall_subgroup)
    except Exception as e:
        print('！！！[Error] fast5 file: ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        signal = or_raw_signal[int(read_start_rel_to_raw):]
        signal_list, signal_mean, signal_std, shift, scale = signal_segmentation(signal,
                                                                                 starts_rel_to_read,
                                                                                 int(event_length[-1]))
        signal_mean = signal_mean / shift
        signal_std = signal_std / scale
    except Exception as e:
        print('[！！！Error] signal: ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        read_bases = np.array(pd.Series(event_bases).apply(get_base_color)) / 300.0
        clean_event_length = np.array(event_length) / 10.0
        x = np.vstack([read_bases, signal_mean, signal_std, clean_event_length, raw_mean, raw_std])
        x = x.T
    except Exception as e:
        print('[！！！Error] input features: ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        y_read, y_qul = get_base_l(default_path, fast5_fn, basecall_tmp_dir, 0, 0, 0)
        out_fasta_fn = args.output_dir + fast5_fn_sg.split('.')[0] + '_out.fasta'
        if not os.path.exists(args.output_dir):
            os.makedirs(args.output_dir)
        prep_read_fasta(fast5_fn, out_fasta_fn, list(y_read))
        shutil.rmtree(basecall_tmp_dir)
        print('[p:::] ' + fast5_fn_sg.split('.')[0] + '_out.fasta was saved......')
    except Exception as e:
        try:
            out_fasta_fn = args.output_dir + fast5_fn_sg.split('.')[0] + '_out.fasta'
            if not os.path.exists(args.output_dir):
                os.makedirs(args.output_dir)
            prep_read_fasta(fast5_fn, out_fasta_fn, list(event_bases))
        except Exception as e:
            print('[！！！Error] stroring : ' + fast5_fn_sg.split('.')[0])


if __name__ == '__main__':
    ar_args=get_args()
    try:
        shutil.rmtree(ar_args.temp_dir)
    except Exception:
        a=1
    check_path(ar_args.temp_dir)
    check_path(ar_args.output_dir)

    fast5_fns = os.listdir(ar_args.fast5_base_dir)
    if int(ar_args.thread) >=len(fast5_fns):
        pool_size = len(fast5_fns)
    else:
        pool_size = int(ar_args.thread)


    # 存在
    # run_time = int(len(fast5_fns)/pool_size)
    run_time = int(len(fast5_fns)/pool_size)
    start_time = time.time()
    p = Pool(pool_size)
    for j in range(run_time):
        for i in range(pool_size):
            # print(pool_size*j+i)
            # print(fast5_fns[pool_size*j+i])
            p.apply_async(provide_fasta,  args=(i, fast5_fns[pool_size*j+i],ar_args,))
    print('[p:::] Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('[s:::] All subprocesses done.')
    # for fast5_fn_sg in fast5_fns:
    #     provide_fasta(fast5_fn_sg, fast5_fn_sg, ar_args,genome_index)
    end_time = time.time()
    print('[s:::] NanoReviser time consuming:%.2f seconds' % (end_time - start_time))

    try:
        shutil.rmtree(ar_args.temp_dir)
    except Exception as e:
        print('！！！[Error] remove tmp dir ' + ar_args.temp_dir + e)


