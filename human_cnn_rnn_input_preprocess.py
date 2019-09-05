# -*- coding: utf-8 -*-
"""
 @File: correcter - cnn_rnn_input_preprocess
 
 @Time: 2019/4/24 2:53 PM
 
 @Author: lotuswang
 
 
"""

import os
from time import sleep
import io
from distutils.version import LooseVersion
from subprocess import call
import re
import h5py
import csv
from itertools import repeat

import numpy as np
import json

import sys
from glob import glob
from collections import defaultdict, namedtuple
import re

from subprocess import call, STDOUT
from tempfile import NamedTemporaryFile
import pandas as pd


class Parser(object):
    def __init__(self):
        self.VERBOSE = True
        self.fast5_base_dir = '../NA12878/test_fast5/'
        self.temp_dir = './tmp/'
        self.train_input_dir = './input/'
        self.temp_file = 'temp_50_reads.fasta'
        self.genome_fn = '../NA12878/NA12878_chr20.fasta'
        self.graphmap_exe = 'graphmap'
        # self.graphmap_exe = '/Users/lotuswang/graphmap/bin/Linux-x64/graphmap'
        # self.bwa_exe = '/Users/lotuswang/bwa/bwa'
        self.failed_reads_filename = 'failed_reads.txt'
        self.basecall_group = 'Basecall_1D_000'
        self.basecall_subgroup = 'BaseCalled_template'
        self.output_format = 'sam'


def parse_fasta(fasta_fn):
    """
    :param fasta_fn:
    :return: fasta_records

    fast_fn:str, the path to Ref.fasta
    fasta_records: dic(),key='Chr', value='tSeq'
    """

    fasta_fp = io.open(fasta_fn, 'r')
    fasta_records = {}
    curr_id = None
    curr_seq = ''
    for line in fasta_fp.readlines():
        if line.startswith('>'):
            if curr_id is not None and curr_seq is not '':
                fasta_records[curr_id] = curr_seq
            curr_seq = ''
            curr_id = line.replace(">", "").strip().split()[0]
        else:
            curr_seq += line.strip()

    # add last record
    if (curr_id is not None and
            curr_seq is not ''):
        fasta_records[curr_id] = curr_seq

    fasta_fp.close()

    return fasta_records


def get_read_data(fast5_fn, basecall_group, basecall_subgroup):
    """
    :param fast5_fn: example:{str} ./input/fast5/id_98490_ch139_read1203_strand.fast5
    :param basecall_group: {str} Default = 'Basecall_1D_000'
    :param basecall_subgroup: {str} Default = 'BaseCalled_template'

    :return:
    **********
    if events.move != 0,save
    if events.move == 2, save 2 elements
    **********
    abs_event_start: event[0]
    start: start times of the bases
    length: bases duration time
    bases:event.model_state[2](move ==1),move == 2
    signal:Raw_signal[abs_event_start:start[-1]+length[-1]]

    all these returns need to be truncated later according to the sam file
    """
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
    except Exception:
        raise NotImplementedError('Error opening file. Likely a corrupted file.')

    try:
        # get albacore version, or if not specified set to 0.0
        albacore_version = LooseVersion(fast5_data['/Analyses/' + basecall_group].attrs['version']
                                        if 'version' in fast5_data['/Analyses/' + basecall_group].attrs
                                        else "0.0")
        if albacore_version <= LooseVersion('0.0'):
            raw_attrs = dict(list(fast5_data['/Raw/Reads/'].values())[0].attrs.items())

            called_dat = fast5_data['/Analyses/' + basecall_group + '/' + basecall_subgroup + '/Events'].value
            called_dat['start'] = called_dat['start'] * 4000 - raw_attrs['start_time']
            called_dat['length'] = called_dat['length'] * 4000
        else:
            called_dat = fast5_data['/Analyses/' + basecall_group + '/' + basecall_subgroup + '/Events'].value
    except Exception:
        raise RuntimeError('No events or corrupted events in file. ' +
                           'Likely a segmentation error .')

    try:
        fastq = fast5_data['/Analyses/' + args.basecall_group + '/' + args.basecall_subgroup + '/Fastq'].value
        quality = str(fastq).split('\\n')
        quality = quality[3]
        quality = re.sub('[\']', '', quality)
        quality = list(quality[2:-2])
        quality = list(map(ord, list(quality)))
        # assert len(quality) == len(called_dat)
    except Exception:
        raise RuntimeError('No fastq file. Likely a basecalling error .')

    start = list()
    bases = list()
    ab_p_model_states = list()
    ab_weights = list()
    for _, start_l, _, _, model_state, move, p_state, weight in zip(called_dat['mean'][::-1],
                                                                    called_dat['start'][::-1],
                                                                    called_dat['stdv'][::-1],
                                                                    called_dat['length'][::-1],
                                                                    called_dat['model_state'][::-1],
                                                                    called_dat['move'][::-1],
                                                                    called_dat['p_model_state'][::-1],
                                                                    called_dat['weights'][::-1]):
        model_state = str(model_state)[2:-1]
        start_l = int(start_l)
        if move == 0:
            continue
        elif move == 1:
            start.append(start_l)
            bases.append(model_state[2])
            ab_p_model_states.append(p_state)
            ab_weights.append(weight)
        elif move == 2:
            start.append(start_l + 2)
            bases.append(model_state[2])
            ab_p_model_states.append(p_state)
            ab_weights.append(weight)
            start.append(start_l)
            bases.append(model_state[1])
            ab_p_model_states.append(p_state)
            ab_weights.append(weight)
        else:
            start.append(start_l)
            bases.append(model_state[2])
            ab_p_model_states.append(p_state)
            ab_weights.append(weight)
    start = start[::-1]
    bases = bases[::-1]
    ab_p_model_states = ab_p_model_states[::-1]
    ab_weights = ab_weights[::-1]

    try:
        length = np.diff(start)
        length = list(length)
        if start[-1] - start[-2] < 5:
            length.append(3.)
        else:
            length.append(5.)
    except Exception:
        raise RuntimeError('Events is too short or ' +
                           'there are too much zero moves.')

    try:
        read_name = list(fast5_data['/Raw/Reads/'].items())[0][0]
        signal = fast5_data[str('/Raw/Reads/') + str(read_name) + '/Signal/'].value
    except Exception:
        raise RuntimeError('No signal stored in the file')

    # TODO: channel_info，raw_attrs are necessary or not？
    # fast5_info = fast5_data['UniqueGlobalKey/channel_id'].attrs
    # raw_attrs = dict(fast5_data['/Raw/Reads/'].values()[0].attrs.items())

    fast5_data.close()
    if len(signal) < int(start[-1] + length[-1]):
        raise RuntimeError('Signal is shorter than the Events')
    else:
        abs_event_start = start[0]
        start = np.array(start) - abs_event_start
        length = np.array(length)
        if len(quality)!=len(ab_p_model_states):
            print(fast5_fn,'quality error')
        return abs_event_start, start, length, bases, quality, signal, ab_p_model_states, ab_weights


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


def prep_graphmap_options(genome_fn, read_fn, out_fn, output_format, num_align_ps):
    """

    :param genome_fn: full path to ref.fasta
    :param read_fn: full path to read.fasta
    :param out_fn: full path to read.sam
    :param output_format: 'sam'
    :param num_align_ps: 1
    :return: list() graphmap_options:  ['align', '-r', genome_fn, '-d', read.fast, '-o', read.sam, '-t', '1']
    """
    return ['align', '-r', genome_fn, '-d', read_fn, '-o', out_fn,
            '-t', str(num_align_ps)]


def parse_fasta(fasta_fn):
    """
    :param fasta_fn:
    :return: fasta_records

    fast_fn:str, the path to Ref.fasta
    fasta_records: dic(),key='Chr', value='tSeq'
    """

    fasta_fp = io.open(fasta_fn, 'r')
    fasta_records = {}
    curr_id = None
    curr_seq = ''
    for line in fasta_fp.readlines():
        if line.startswith('>'):
            if curr_id is not None and curr_seq is not '':
                fasta_records[curr_id] = curr_seq
            curr_seq = ''
            curr_id = line.replace(">", "").strip().split()[0]
        else:
            curr_seq += line.strip()

    # add last record
    if (curr_id is not None and
            curr_seq is not ''):
        fasta_records[curr_id] = curr_seq

    fasta_fp.close()
    return fasta_records


def align_to_genome(out_fn, graphmap_exe, mapper_options):
    """

    :param graphmap_exe: full path to graphmap，should like ~/graphmap/bin/Linux-x64/graphmap
    :param mapper_options: ['align', '-r', genome_fn, '-d', read.fast, '-o', read.sam, '-t', '1']
    :return: dict() sam_record: key = ('qName', 'flag', 'rName', 'pos', 'mapq',
                                        'cigar', 'rNext', 'pNext', 'tLen', 'seq', 'qual')
    """
    FNULL = open(os.devnull, 'w')
    stdout_sink = FNULL
    exitStatus = call([graphmap_exe, ] + mapper_options, stdout=stdout_sink, stderr=FNULL)
    FNULL.close()
    if exitStatus == 0:
        with open(out_fn, 'r') as out_fp:
            out_fp.seek(0)
            align_output = out_fp.readlines()
    else:
        raise RuntimeError('Align Error, please check your graphmap or bwa mem')

    SAM_FIELDS = ('qName', 'flag', 'rName', 'pos', 'mapq',
                  'cigar', 'rNext', 'pNext', 'tLen', 'seq', 'qual')

    r_sam_record = dict()
    for line in align_output:
        if line.startswith('@'):
            continue
        else:
            r_sam_record = dict(zip(SAM_FIELDS, line.strip().split()))
    if not r_sam_record:
        raise RuntimeError('Map Error, there is no read record in the sam file')
    elif len(r_sam_record) < len(SAM_FIELDS) or r_sam_record['rName'] == '*':
        raise RuntimeError('Map Error, the read is unmapped.')
    else:
        return r_sam_record


def comp_base(base):
    # replace non-ACGT bases with dash
    COMP_BASES = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-', 'N': 'N'}
    try:
        return COMP_BASES[base]
    except KeyError:
        return 'N'


def rev_comp(seq):
    return ''.join(comp_base(b) for b in seq[::-1])


def parse_sam_record(r_sam_record, genome_index):
    # parse cigar string
    CIGAR_PAT = re.compile('(\d+)([MIDNSHP=X])')
    genomeLoc = dict()
    cigar = [
        (int(reg_len), reg_type) for reg_len, reg_type in
        CIGAR_PAT.findall(r_sam_record['cigar'])]
    if len(cigar) < 1:
        raise RuntimeError('Invalid cigar string produced.')

    strand = '-' if int(r_sam_record['flag']) & 0x10 else '+'
    if strand == '-':
        cigar = cigar[::-1]
    # cigar_or = cigar
    # record clipped bases and remove from query seq as well as cigar
    qSeq = r_sam_record['seq'] if strand == '+' else rev_comp(
        r_sam_record['seq'])
    qSeq = str(qSeq)
    start_clipped_bases = 0
    end_clipped_bases = 0
    # handle clipping elements (H and S)
    if cigar[0][1] == 'H':
        start_clipped_bases += cigar[0][0]
        cigar = cigar[1:]
    if cigar[-1][1] == 'H':
        end_clipped_bases += cigar[-1][0]
        cigar = cigar[:-1]
    if cigar[0][1] == 'S':
        start_clipped_bases += cigar[0][0]
        qSeq = qSeq[cigar[0][0]:]
        cigar = cigar[1:]
    if cigar[-1][1] == 'S':
        end_clipped_bases += cigar[-1][0]
        qSeq = qSeq[:-cigar[-1][0]]
        cigar = cigar[:-1]

    tLen = sum([reg_len for reg_len, reg_type in cigar
                if reg_type in 'MDN=X'])
    tSeq = genome_index[r_sam_record['rName']][
        int(r_sam_record['pos']) - 1:
        int(r_sam_record['pos']) + tLen - 1]
    tSeq = str(tSeq)
    if strand == '-':
        tSeq = rev_comp(tSeq)

    # check that cigar starts and ends with matched bases
    while cigar[0][1] not in 'M=X':
        if cigar[0][1] in 'IP':
            tSeq = tSeq[cigar[0][0]:]
        else:
            qSeq = qSeq[cigar[0][0]:]
            start_clipped_bases += cigar[0][0]
        cigar = cigar[1:]
    while cigar[-1][1] not in 'M=X':
        if cigar[-1][1] in 'IP':
            tSeq = tSeq[:-cigar[-1][0]]
        else:
            qSeq = qSeq[:-cigar[-1][0]]
            end_clipped_bases += cigar[0][0]
        cigar = cigar[:-1]

    qLen = sum([reg_len for reg_len, reg_type in cigar
                if reg_type in 'MIP=X'])
    assert len(qSeq) == qLen , 'Read sequence from SAM and cooresponding cigar string do not agree.'

    # create pairwise alignment via zipped pairs
    readVals = []
    refVals = []
    mapVals = []
    for reg_len, reg_type in cigar:
        if reg_type in 'M=X':
            tmp_read = list(qSeq[:reg_len])
            readVals.extend(tmp_read)
            tmp_ref = list(tSeq[:reg_len])
            refVals.extend(tmp_ref)
            for read, ref in zip(tmp_read, tmp_ref):
                if read == ref:
                    mapVals.extend('M')
                else:
                    mapVals.extend('X')
            qSeq = qSeq[reg_len:]
            tSeq = tSeq[reg_len:]
        elif reg_type in 'IP':
            tmp_read = list(qSeq[:reg_len])
            readVals.extend(tmp_read)
            tmp_ref = list('-'*reg_len)
            refVals.extend(tmp_ref)
            mapVals.extend(list('I'*reg_len))
            qSeq = qSeq[reg_len:]
        else:
            tmp_ref = list(tSeq[:reg_len])
            refVals.extend(tmp_ref)
            tmp_read = list('-'*reg_len)
            readVals.extend(tmp_read)
            mapVals.extend(list('D' * reg_len))
            tSeq = tSeq[reg_len:]

    genomeLoc['Start'] = int(r_sam_record['pos']) - 1
    genomeLoc['Strand'] = strand
    genomeLoc['Chrom'] = r_sam_record['rName']
    return readVals, refVals, mapVals, genomeLoc, start_clipped_bases, end_clipped_bases


def fix_raw_starts_for_clipped_bases(start_clipped_bases, end_clipped_bases,
                                     starts_rel_to_read, event_length, read_start_rel_to_raw,
                                     ab_p_model_states, ab_weights, quality):
    if start_clipped_bases > 0:
        # print(start_clipped_bases)
        start_clipped_bases = int(start_clipped_bases)
        start_clipped_obs = int(starts_rel_to_read[int(start_clipped_bases)])
        ab_p_model_states = ab_p_model_states[int(start_clipped_bases):]
        ab_weights = ab_weights[int(start_clipped_bases):]
        quality = quality[int(start_clipped_bases):]
        event_length = event_length[int(start_clipped_bases):]
        starts_rel_to_read = starts_rel_to_read[start_clipped_bases:] - start_clipped_obs
        read_start_rel_to_raw += start_clipped_obs
        read_start_rel_to_raw = int(read_start_rel_to_raw)

    if end_clipped_bases > 0:
        # print(end_clipped_bases)
        end_clipped_bases = int(end_clipped_bases)
        starts_rel_to_read = starts_rel_to_read[:-1 * end_clipped_bases]
        ab_p_model_states = ab_p_model_states[:-1 * end_clipped_bases]
        ab_weights = ab_weights[:-1 * end_clipped_bases]
        quality = quality[:-1 * end_clipped_bases]
        event_length = event_length[:-1 * end_clipped_bases]

    return starts_rel_to_read, event_length, read_start_rel_to_raw, ab_p_model_states, ab_weights,quality


def clean_read_map_ref(readVals, mapVals, refVals):
    clean_readVals = list()
    clean_mapVals = list()
    clean_refVals = list()
    del_num_one = 0
    del_num_more = 0
    for read_v1, map_v1, ref_v1, read_v2, map_v2, ref_v2 in \
            zip(readVals[:-1], mapVals[:-1], refVals[:-1], readVals[1:], \
                mapVals[1:], refVals[1:]):
        if map_v1 in 'MXI' and map_v2 in 'MXI':
            clean_readVals.append(read_v1)
            clean_mapVals.append(map_v1)
            clean_refVals.append(ref_v1)
            continue
        elif map_v1 in 'MXI' and map_v2 == 'D':
            clean_readVals.append(read_v1)
            clean_mapVals.append(map_v2)
            clean_refVals.append('D')
            del_num_one += 1
            continue
        elif map_v1 == 'D' and map_v2 == 'D':
            del_num_more += 1
            continue
        elif map_v1 == 'D' and map_v2 in 'MXI':
            continue
        else:
            print(read_v1, map_v1, ref_v1, read_v2, map_v2, ref_v2)
    clean_readVals.append(readVals[-1])
    clean_mapVals.append(mapVals[-1])
    clean_refVals.append(refVals[-1])
    # print(del_num_one, del_num_more)
    assert len(clean_readVals) == len(clean_refVals)
    assert len(clean_readVals) == len(clean_mapVals)
    return clean_readVals, clean_mapVals, clean_refVals


def signal_segmentation(signal, starts_rel_to_read, last_dur, query_len=50):
    query_len = int(query_len)
    if query_len % 2 == 0:
        query_ahead = query_len / 2
        query_tail = query_len / 2
    else:
        query_len = (query_len - 1)
        query_ahead = query_len / 2
        query_tail = 1 + query_len / 2

    tmp_raw_signal = np.array(signal).astype(float)
    shift = np.median(tmp_raw_signal)
    scale = np.median(np.abs(tmp_raw_signal - shift))
    signal = (raw_signal - shift) / (scale+0.0001)
    signal_mean = list()
    signal_std = list()
    st_e = int(starts_rel_to_read[-1])
    en_e = st_e + int(last_dur)
    # print(st_e, en_e)
    signal_list = list()
    signal_len = len(tmp_raw_signal)
    try:
        for st, en in zip(starts_rel_to_read[:-1], starts_rel_to_read[1:]):
            st = int(st)
            en = int(en)
            tmp_sig = signal[st:en]
            if st - query_ahead <= 0:
                tmp_st = 0
            else:
                tmp_st = int(st - query_ahead)
            if st + query_tail >= signal_len:
                tmp_en = signal_len
            else:
                tmp_en = int(st + query_tail)
            tmp1_sig = raw_signal[tmp_st:tmp_en]
            pad_len = query_len - len(tmp1_sig)
            if pad_len > 0:
                if pad_len % 2:
                    pad_len = int(pad_len / 2)
                    tmp1_sig = list(np.pad(tmp1_sig, (pad_len + 1, pad_len),
                                           'constant',
                                           constant_values=(0, 0)))
                else:
                    pad_len = int((pad_len + 1) / 2)
                    tmp1_sig = list(np.pad(tmp1_sig, (pad_len, pad_len),
                                           'constant',
                                           constant_values=(0, 0)))
            tmp1_sig = np.array(tmp1_sig).astype(float) / (shift+0.0001)
            signal_list.append(tmp1_sig)
            tmp_mean = np.mean(tmp_sig).astype(float)
            tmp_std = np.std(tmp_sig).astype(float)
            signal_mean.append(tmp_mean)
            signal_std.append(tmp_std)

        tmp_sig = signal[st_e:en_e]
        if st_e - query_ahead <= 0:
            tmp_st = 0
        else:
            tmp_st = int(st_e - query_ahead)
        if st_e + query_tail >= signal_len:
            tmp_en = int(signal_len)
        else:
            tmp_en = int(st_e + query_tail)
        tmp1_sig = raw_signal[tmp_st:tmp_en]
        pad_len = query_len - len(tmp1_sig)
        if pad_len > 0:
            if pad_len % 2:
                pad_len = int(pad_len / 2)
                tmp1_sig = list(np.pad(tmp1_sig, (pad_len + 1, pad_len),
                                       'constant',
                                       constant_values=(0, 0)))
            else:
                pad_len = int((pad_len + 1) / 2)
                tmp1_sig = list(np.pad(tmp1_sig, (pad_len, pad_len),
                                       'constant',
                                       constant_values=(0, 0)))
        tmp1_sig = np.array(tmp1_sig).astype(float) / scale
        signal_list.append(tmp1_sig)
        tmp_mean = np.mean(tmp_sig)
        tmp_std = np.std(tmp_sig)
        signal_mean.append(tmp_mean)
        signal_std.append(tmp_std)
        assert len(signal_mean) == len(signal_std)
    except Exception:
        raise RuntimeError('Signal segmentation Error')
    return np.array(signal_list), signal_mean, signal_std, shift, scale


def get_base_color(base):
    base_to_color = {'A':250, 'G':180, 'T':100, 'C':30}
    return base_to_color.get(base,0)


def get_base_label(base):
    base_to_label = {'A':5, 'G':4, 'T':3, 'C':2,'-':1, 'D':0}
    return base_to_label.get(base,0)


if __name__ == '__main__':
    args = Parser()
    genome_fn = args.genome_fn
    genome_index = parse_fasta(genome_fn)
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
            print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.fasta was saved for mapping......')

            num_align_ps = 1
            out_fn = args.temp_dir + fast5_fn_sg.split('.')[0] + '.sam'
            graphmap_options = prep_graphmap_options(args.genome_fn, read_fasta_fn, out_fn, args.output_format,
                                                 num_align_ps)


            sam_records = align_to_genome(out_fn, args.graphmap_exe, graphmap_options)
            print('[p:::] ' + fast5_fn_sg.split('.')[0] + '.sam has been loaded......')

        except Exception as e:
            print('！！！[Error] ' + fast5_fn_sg.split('.')[0] + str(e))
            continue
        try:
            readVals, refVals, mapVals, genomeLoc, \
            start_clipped_bases, end_clipped_bases = parse_sam_record(sam_records,
                                                                      genome_index)
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
            print('[s:::] ' + fast5_fn_sg.split('.')[0] + '.npz has been saved......')
            os.remove(out_fn)
            os.remove(read_fasta_fn)
            os.remove('../train_human/NA12878_chr20.fasta.gmidx')

        except Exception as e:
            print('[！！！Error] ' + fast5_fn_sg.split('.')[0] + str(e))
            continue














