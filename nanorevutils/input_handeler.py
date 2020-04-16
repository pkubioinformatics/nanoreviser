# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - input_handeler
 
 @Time: 2019/4/14 4:52 PM
 
 @Author: lotuswang
 
 
"""
import io
import re


def comp_base(base):
    # replace non-ACGT bases with dash
    COMP_BASES = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '-': '-', 'N': 'N'}
    try:
        return COMP_BASES[base]
    except KeyError:
        return 'N'


def rev_comp(seq):
    return ''.join(comp_base(b) for b in seq[::-1])


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


