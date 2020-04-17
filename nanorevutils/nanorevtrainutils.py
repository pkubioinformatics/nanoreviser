# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - nanorevtrainutils
 
 @Time: 2020/4/16 12:53 PM
 
 @Author: lotuswang
 
 
"""

import warnings

warnings.filterwarnings('ignore')

import os
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
from nanorevutils.alignutils import align_to_genome, parse_sam_record, prep_graphmap_options
from nanorevutils.preprocessing import clean_read_map_ref, fix_raw_starts_for_clipped_bases

NB_CLASS = 6
SENT_LEN = 13
SIGNEL_LEN = 50
VEC_LEN = 6

def prep_read_fasta(fast5_fn, read_fasta_fn, bases):
    """
    :param fast5_fn: ./input/fast5/id_98490_ch139_read1203_strand.fast5
    :param read_fasta_fn: ./input/tmp/id_98490_ch139_read1203_strand.fasta
    :param bases: ['G', 'T', 'T', 'G', 'C', 'T', 'T', 'C', 'G', 'T', 'T']
    :return:
    return None
    """
    try:
        # reads_fasta = ''
        reads_fasta = ">" + fast5_fn.replace(' ', '|||') + '\n' + ''.join(bases) + '\n'
        # print(reads_fasta)
        read_fp = open(read_fasta_fn, 'w')
        read_fp.write(reads_fasta)
        read_fp.close()
    except Exception as e:
        raise NotImplementedError('Error in writing .fasta file')
    return True


def handel_input_fast5(fast5_fn_sg, args, genome_index):
    fast5_fn = os.path.join(args.fast5_base_dir, fast5_fn_sg)
    try:
        (read_start_rel_to_raw, starts_rel_to_read, event_length, event_bases,
         or_raw_signal, raw_mean, raw_std) = get_read_data(fast5_fn, args.basecall_group, args.basecall_subgroup)
    except Exception as e:
        raise EOFError('！！！[Error] ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        read_fasta_fn = args.temp_dir + fast5_fn_sg.split('.')[0] + '.fasta'
        prep_read_fasta(fast5_fn, read_fasta_fn, event_bases)
        print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.fasta was saved for mapping......')

        num_align_ps = 1
        out_fn = args.temp_dir + fast5_fn_sg.split('.')[0] + '.sam'
        graphmap_options = prep_graphmap_options(args.genome_fn, read_fasta_fn, out_fn, args.output_format,
                                                 num_align_ps)

        sam_records = align_to_genome(out_fn, args.graphmap_exe, graphmap_options)
        print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.sam has been loaded......')

    except Exception as e:
        raise RuntimeError('！！！[Error] ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        readVals, refVals, mapVals, genomeLoc, \
        start_clipped_bases, end_clipped_bases = parse_sam_record(sam_records,
                                                                  genome_index)
        print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.sam is mapping......')

        starts_rel_to_read, event_length, read_start_rel_to_raw, raw_mean, raw_std = \
            fix_raw_starts_for_clipped_bases(int(start_clipped_bases),
                                             int(end_clipped_bases),
                                             starts_rel_to_read,
                                             event_length,
                                             int(read_start_rel_to_raw),
                                             raw_mean,
                                             raw_std,
                                             )
        clean_readVals, clean_mapVals, clean_refVals, clean_refVals2, \
            = clean_read_map_ref(readVals,
                                 mapVals,
                                 refVals,)
        signal = or_raw_signal[int(read_start_rel_to_raw):]
        signal_list, signal_mean, signal_std, shift, scale = \
            signal_segmentation(signal, starts_rel_to_read, int(event_length[-1]))

        readvals = np.array(pd.Series(clean_readVals).apply(get_base_color))
        refvals = np.array(pd.Series(clean_refVals).apply(get_base_label))
        refvals2 = np.array(pd.Series(clean_refVals2).apply(get_base_label))

        save_name = str(args.train_input_dir) + str(fast5_fn_sg).split('.')[0]
        starts = np.array(starts_rel_to_read)
        signal_len = np.array(event_length)
        starts = np.array(starts_rel_to_read)
        signal_len = np.array(event_length)
        np.savez(save_name,
                 refvals=np.array(refvals),
                 refvals2=np.array(refvals2),
                 readVals=np.array(readvals),
                 signal_mean=signal_mean,
                 signal_std=signal_std,
                 signal_len=np.array(signal_len),
                 ab_mean=np.array(raw_mean),
                 ab_std=np.array(raw_std),
                 signal_x=signal_list,
                 mapvals=np.array(clean_mapVals),
                 starts=np.array(starts),
                 scale=scale,
                 shift=shift)
        if not args.test_mode:
            print('[s:::] ' + fast5_fn_sg.split('.')[0] + '.npz has been saved......')
        os.remove(out_fn)
        os.remove(read_fasta_fn)

    except Exception as e:
        raise RuntimeError('[！！！Error] ' + fast5_fn_sg.split('.')[0] + str(e))


def train_preprocessing(args):
    genome_fn = args.genome_fn
    genome_index = parse_fasta(genome_fn)
    if not args.test_mode:
        print(genome_fn, 'has been load......')
    if args.read_counts==0 or args.read_counts>=len(os.listdir(args.fast5_base_dir)):
        fast5_fns = os.listdir(args.fast5_base_dir)
    else:
        fast5_fns = os.listdir(args.fast5_base_dir)[:int(args.read_counts)]
    for fast5_fn_sg in fast5_fns:
        handel_input_fast5(fast5_fn_sg, args, genome_index)


def get_trainning_input(test_mode, train_input_dir, window_size=13):
    train_list = os.listdir(train_input_dir)
    x=[0]
    i=0
    SENT_LEN = window_size
    for train_fn in train_list:
        try:
            train_fn = os.path.join(train_input_dir, train_fn)
            if train_fn.endswith('.npz'):
                npzfile = np.load(train_fn)
                shift = npzfile['shift']
                scale = npzfile['scale']
                readvals = npzfile['readVals']
                signal_mean = npzfile['signal_mean'] / shift
                signal_std = npzfile['signal_std'] / scale
                signal_len = npzfile['signal_len'] / 10.0
                ab_mean = npzfile['ab_mean']
                ab_std = npzfile['ab_std']
                # tmp_read = np.array(pd.Series(readvals).apply(get_base_color)) / 300.0
                tmp_read = readvals/300.0
                tmp_x = np.vstack([tmp_read, signal_mean, signal_std, signal_len, ab_mean, ab_std])
                tmp_signal_x = npzfile['signal_x']
                tmp_y = npzfile['refvals']
                tmp_y2 = npzfile['refvals2']
                i += 1
                if i < 2:
                    read_x = readvals
                    x = tmp_x
                    signal_x = tmp_signal_x
                    y = tmp_y
                    y2 = tmp_y2
                else:
                    x = np.append(x, tmp_x, axis=-1)
                    signal_x = np.append(signal_x, tmp_signal_x, axis=0)
                    y = np.append(y, tmp_y, axis=0)
                    y2 = np.append(y2, tmp_y2, axis=0)
                    read_x = np.append(read_x, readvals, axis=0)
        except Exception as e:
            print('！！！[Error] training input file:', train_fn)
            continue
    x = x.T
    assert len(x)>2*SENT_LEN
    try:
        x = np.array(x, dtype='float')
        # y = np.array(pd.Series(y).apply(get_base_label))
        y = np.array(y, dtype='float')
        # y2 = np.array(pd.Series(y2).apply(get_base_label))
        y2 = np.array(y2, dtype='float')
        x_train_tmp = list()
        for i in range(len(x) - SENT_LEN):
            tmp_x = x[i:i + SENT_LEN]
            x_train_tmp.append(tmp_x)
        x_train = np.array(x_train_tmp)
        del x_train_tmp
        signal_x_train_tmp = list()
        for i in range(len(signal_x) - SENT_LEN):
            tmp_x = signal_x[i:i + SENT_LEN]
            signal_x_train_tmp.append(tmp_x)
        # x_train = np.array(tmp_x)p
        # print(len(signal_x_train_tmp))
        signal_x_train = np.array(signal_x_train_tmp, dtype='float')
        SET_BEF = int((SENT_LEN - 1) / 2)
        SET_AFT = int((SENT_LEN + 1) / 2)
        y_train = np.reshape(y[SET_BEF:-SET_AFT], [-1, 1, ])
        y_train2 = np.reshape(y2[SET_BEF:-SET_AFT] - 1, [-1, 1, ])
        if not test_mode:
            print('[p:::] input files has been load......')
        return x_train, signal_x_train, y_train, y_train2
    except Exception as e:
        raise RuntimeError('！！！[Error] fatal errors in loading training data.')


