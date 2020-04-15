# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - output_handeler
 
 @Time: 2019/4/14 4:46 PM
 
 @Author: lotuswang
 
 
"""
import os
from nanorevutils.nanorev_fast5_handeler import handle_opts, opt_main
from nanorevutils.lstmmodel import get_model1, get_model2


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
                dna_seq=fastq_output[1]
                dna_qul=fastq_output[3]
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




