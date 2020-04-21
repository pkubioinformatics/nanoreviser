# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - check_nanoreviser
 
 @Time: 2020/4/17 2:56 PM
 
 @Author: lotuswang
 
 
"""

import os
import shutil

def delet_files(log_fn):
    try:
        os.remove(log_fn)
    except Exception as e:
        a=1
    if os.path.exists('../unitest_training_results/'):
        shutil.rmtree('../unitest_training_results')

if __name__ == '__main__':
    log_fn='./unitest/unitest_log.txt'
    tag = 0
    with open(log_fn, 'r') as lg:
        result_list = lg.readlines()
    for result in result_list:
        tags=result.split(' - ')
        tag1 = tags[2]
        tag2 = tags[-1].split(',')[0]
        if not (tag1=='INFO' and tag2=='Congratulations'):
            tag = 1
            delet_files(log_fn)
            print('Error!!! NanoReviser is not properly installed ! plese check your platform! ')
            break

    if tag==0:
        delet_files(log_fn)
        print("Congratulations, please have fun with NanoReviser :)")

