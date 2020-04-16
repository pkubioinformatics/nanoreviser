# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - preprocessing
 
 @Time: 2019/10/7 10:41 AM
 
 @Author: lotuswang
 
 
"""

import warnings
warnings.filterwarnings('ignore')
import numpy as np



def fix_raw_starts_for_clipped_bases(start_clipped_bases, end_clipped_bases,
                                     starts_rel_to_read, event_length, read_start_rel_to_raw,
                                     ab_p_model_states, ab_weights):
    if start_clipped_bases > 0:
        # print(start_clipped_bases)
        start_clipped_bases = int(start_clipped_bases)
        start_clipped_obs = int(starts_rel_to_read[int(start_clipped_bases)])
        ab_p_model_states = ab_p_model_states[int(start_clipped_bases):]
        ab_weights = ab_weights[int(start_clipped_bases):]
        # quality = quality[int(start_clipped_bases):]
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
        # quality = quality[:-1 * end_clipped_bases]
        event_length = event_length[:-1 * end_clipped_bases]

    return starts_rel_to_read, event_length, read_start_rel_to_raw, ab_p_model_states, ab_weights


def clean_read_map_ref(readVals, mapVals, refVals):
    clean_readVals = list()
    clean_mapVals = list()
    clean_refVals = list()
    clean_refVals2 = list()
    del_num_one = 0
    del_num_more = 0
    for read_v1, map_v1, ref_v1, read_v2, map_v2, ref_v2 in \
            zip(readVals[:-1], mapVals[:-1], refVals[:-1], readVals[1:], \
                mapVals[1:], refVals[1:]):
        if map_v1 in 'MXI' and map_v2 in 'MXI':
            clean_readVals.append(read_v1)
            clean_mapVals.append(map_v1)
            clean_refVals.append(ref_v1)
            clean_refVals2.append(ref_v1)
            continue
        elif map_v1 in 'MXI' and map_v2 == 'D':
            clean_readVals.append(read_v1)
            clean_mapVals.append(map_v2)
            clean_refVals.append('D')
            clean_refVals2.append(ref_v1)
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
    clean_refVals2.append(refVals[-1])
    # print(del_num_one, del_num_more)
    assert len(clean_readVals) == len(clean_refVals)
    assert len(clean_readVals) == len(clean_mapVals)
    return clean_readVals, clean_mapVals, clean_refVals, clean_refVals2


def signal_segmentation(raw_signal, starts, last_dur, query_len=50):
    query_len = int(query_len)
    if query_len % 2 == 0:
        query_ahead = query_len / 2
        query_tail = query_len / 2
    else:
        query_len = (query_len - 1)
        query_ahead = query_len / 2
        query_tail = 1 + query_len / 2

    tmp_raw_signal = np.array(raw_signal).astype(float)
    shift = np.median(tmp_raw_signal)
    scale = np.median(np.abs(tmp_raw_signal - shift))
    # signal = (raw_signal - shift) / (scale+0.0001)
    signal_mean = list()
    signal_std = list()
    st_e = int(starts[-1])
    en_e = st_e + int(last_dur)
    signal_list = list()
    signal_len = len(raw_signal)
    # print(st_e, en_e)
    try:
        for st, en in zip(starts[:-1], starts[1:]):
            st = int(st)
            en = int(en)
            tmp_sig = raw_signal[st:en]
            if st - query_ahead <= 0:
                tmp_st = 0
            else:
                tmp_st = int(st - query_ahead)
            if st + query_tail >= signal_len:
                tmp_en = signal_len
            else:
                tmp_en = int(st + query_tail)
            tmp1_sig = (np.array(raw_signal[tmp_st:tmp_en]) - shift) / scale
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
            tmp1_sig = np.array(tmp1_sig)
            signal_list.append(tmp1_sig)
            tmp_mean = np.mean(1.0*tmp_sig)
            tmp_std = np.std(1.0*tmp_sig)
            signal_mean.append(tmp_mean)
            signal_std.append(tmp_std)

        tmp_sig = raw_signal[st_e:en_e]
        if st_e - query_ahead <= 0:
            tmp_st = 0
        else:
            tmp_st = int(st_e - query_ahead)
        if st_e + query_tail >= signal_len:
            tmp_en = int(signal_len)
        else:
            tmp_en = int(st_e + query_tail)
        tmp1_sig = (np.array(raw_signal[tmp_st:tmp_en]) - shift) / scale
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
        tmp1_sig = np.array(tmp1_sig)
        signal_list.append(tmp1_sig)
        tmp_mean = np.mean(1.0*tmp_sig)
        tmp_std = np.std(1.0*tmp_sig)
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




