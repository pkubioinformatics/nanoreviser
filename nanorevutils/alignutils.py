# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - alignutils
 
 @Time: 2019/4/14 4:55 PM
 
 @Author: lotuswang
 
 
"""
import os
import re
from subprocess import call, STDOUT


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