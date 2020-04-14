# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - input_handeler
 
 @Time: 2019/4/14 4:52 PM
 
 @Author: lotuswang
 
 
"""
import io

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


