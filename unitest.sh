#########################################################################
# File Name:unitest.sh
# Author:Lotus Wang
# Description: unitest NanoReviser
#########################################################################

#!/bin/bash
#set -x

# 1. Test NanoReviser to generate fasta or fastq
# ---------------------------

python NanoReviser.py -d ./unitest/test_data/fast5/ -o ./fastq_file/ -F fasta --disable_print --test_mode
python NanoReviser.py -d ./unitest/test_data/fast5/ -o ./fastq_file/ -F fastq --disable_print --test_mode

# 1. Test NanoReviser_train for generate model file
# ---------------------------

pyton NanoReviser_train.py -d ./unitest/training_data/fast5/ -r ./unitest/training_data/reference.fasta -o ./unintest/training_results/ -m unitest --disable_print --test_mode
python ./unitest/check_nanoreviser_train.py
