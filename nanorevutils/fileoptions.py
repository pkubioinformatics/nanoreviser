# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - fileoptions
 
 @Time: 2020/4/12 4:58 PM
 
 @Author: lotuswang
 
 
"""
import os
import sys
from pathlib import Path
import shutil
import json


def this_folder():
    """ We need this in every script as the scripts may not live anywhere near
    the albacore package itself.

    """
    if getattr(sys, 'frozen', False):
        # The application is frozen
        return os.path.dirname(sys.executable)
    else:
        # The application is not frozen
        return os.path.dirname(__file__)


def check_path(my_path):
    my_file = Path(str(my_path))
    my_path=str(my_path)
    if not os.path.exists(my_file):
        try:
            os.makedirs(my_path)
        except Exception as e:
            raise FileNotFoundError('！！！[Error] make dir ' + my_path + e)


def copy_file(fn, dstpath='./tmp/basecall_tmp/'):
    src = fn
    if not os.path.exists(dstpath):
        os.makedirs(dstpath)
    if os.path.exists(src) and os.path.exists(dstpath):
        #src = os.path.join(srcpath, file_name)
        result = shutil.copy(src, dstpath)
        if result:
            # print('[running] ',fn, 'has been copied.....')
            a =1
        else:
            print('!!![Error] ',fn, 'cannot move.....')
    else:
        print('!!![Error] ',src, 'does not exit')


def dict_to_json_write_file(hash_set, save_fn):
    with open(save_fn, 'w') as f:
        json.dump(hash_set, f)


def json_file_to_dict(load_fn):
    with open(load_fn, 'r') as f:
        hash_set = json.load(f)
    return hash_set