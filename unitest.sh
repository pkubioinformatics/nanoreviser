#########################################################################
# File Name:unitest.sh
# Author:Lotus Wang
# Description: unitest NanoReviser
#########################################################################

#!/bin/bash
#set -x

# 1. Test NanoReviser to generate fasta or fastq
# ---------------------------

python NanoReviser.py -d ./unitest/test_data/fast5/ -o ./unitest/fasta_file/ -F fasta --test_mode
python NanoReviser.py -d ./unitest/test_data/fast5/ -o ./unitest/fastq_file/ -F fastq --test_mode

# 2. Test NanoReviser_train for generate model file
# ---------------------------

# python NanoReviser_train.py -d ./unitest/training_data/fast5/ -r ./unitest/training_data/reference.fasta -o ./unitest_training_results/ -m unitest --test_mode
python NanoReviser_train.py  --test_mode

# 3. check the output of nanoviser
# ---------------------------
# python ./unitest/check_nanoreviser.py

echo "unitest has been completed."
