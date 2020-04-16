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


def train_preprocessing(args):
    if args.test_mode:
        logger = logger_config(log_path='./unitest/unitest_log.txt', logging_name='unitest')
    genome_fn = args.genome_fn
    genome_index = parse_fasta(genome_fn)
    if not args.test_mode:
        print(genome_fn, 'has been load......')
    fast5_fns = os.listdir(args.fast5_base_dir)
    for fast5_fn_sg in fast5_fns:
        fast5_fn = os.path.join(args.fast5_base_dir, fast5_fn_sg)
        try:
            (read_start_rel_to_raw, starts_rel_to_read, event_length, event_bases, raw_quality,
             raw_signal, ab_p_model_states, ab_weights) = get_read_data(fast5_fn,
                                                                        args.basecall_group,
                                                                        args.basecall_subgroup)
        except Exception as e:
            print('！！！[Error] ' + fast5_fn_sg.split('.')[0] + str(e))
            continue
        try:
            read_fasta_fn = args.temp_dir + fast5_fn_sg.split('.')[0] + '.fasta'
            prep_read_fasta(fast5_fn, read_fasta_fn, event_bases)
            if not args.test_mode:
                print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.fasta was saved for mapping......')

            num_align_ps = 1
            out_fn = args.temp_dir + fast5_fn_sg.split('.')[0] + '.sam'
            graphmap_options = prep_graphmap_options(args.genome_fn, read_fasta_fn, out_fn, args.output_format,
                                                     num_align_ps)

            sam_records = align_to_genome(out_fn, args.graphmap_exe, graphmap_options)
            if not args.test_mode:
                print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.sam has been loaded......')

        except Exception as e:
            print('！！！[Error] ' + fast5_fn_sg.split('.')[0] + str(e))
            continue
        try:
            readVals, refVals, mapVals, genomeLoc, \
            start_clipped_bases, end_clipped_bases = parse_sam_record(sam_records,
                                                                      genome_index)
            if not args.test_mode:
                print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.sam is mapping......')

            starts_rel_to_read, event_length, read_start_rel_to_raw, ab_p_model_states, ab_weights, quality = \
                fix_raw_starts_for_clipped_bases(int(start_clipped_bases),
                                                 int(end_clipped_bases),
                                                 starts_rel_to_read,
                                                 event_length,
                                                 int(read_start_rel_to_raw),
                                                 ab_p_model_states,
                                                 ab_weights,
                                                 raw_quality)
            clean_readVals, clean_mapVals, clean_refVals = clean_read_map_ref(readVals, mapVals, refVals)
            signal = raw_signal[int(read_start_rel_to_raw):]
            signal_list, signal_mean, signal_std, shift, scale = \
                signal_segmentation(signal, starts_rel_to_read, int(event_length[-1]))

            readvals = np.array(pd.Series(clean_readVals).apply(get_base_color))
            refvals = np.array(pd.Series(clean_refVals).apply(get_base_label))

            save_name = str(args.train_input_dir) + str(fast5_fn_sg).split('.')[0]
            starts = np.array(starts_rel_to_read)
            signal_len = np.array(event_length)
            np.savez(save_name,
                     readvals=readvals,
                     clean_readVals=np.array(clean_readVals),
                     signal_mean=signal_mean,
                     signal_std=signal_std,
                     signal_len=signal_len,
                     ab_p_model_states=np.array(ab_p_model_states),
                     ab_weights=np.array(ab_weights),
                     signal_x=signal_list,
                     refvals=refvals,
                     mapvals=np.array(clean_mapVals),
                     starts=starts,
                     quality=quality)
            if not args.test_mode:
                print('[s:::] ' + fast5_fn_sg.split('.')[0] + '.npz has been saved......')
            os.remove(out_fn)
            os.remove(read_fasta_fn)

        except Exception as e:
            print('[！！！Error] ' + fast5_fn_sg.split('.')[0] + str(e))
            continue