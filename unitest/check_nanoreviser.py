# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - check_nanoreviser
 
 @Time: 2020/4/17 2:56 PM
 
 @Author: lotuswang
 
 
"""

import os
import shutil

if __name__ == '__main__':
    log_fn='./unitest/unitest_log.txt'
    with open(log_fn, 'r') as lg:
        result_list = lg.readlines()
    assert len(result_list)==3
    for result in result_list:
        tags=result.split(' - ')
        tag1 = tags[2]
        tag2 = tags[-1].split(',')[0]
        assert (tag1=='INFO' and tag2=='Congratulations')
    try:
        os.remove(log_fn)
    except Exception as e:
        a=1
    if os.path.exists('../unitest_training_results/'):
        shutil.rmtree('../unitest_training_results')
    print("Congratulations, please have fun with NanoReviser :)")

