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
import time
from albacore.path_utils import get_default_path


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

def model_fn_generate(args, model_tag):
    model_train_fn = args.train_model_dir + \
                 'train_'+str(args.species) + \
                 '_win' + str(args.window_size) + '_' + \
                 str(args.epochs) +'ep_' + \
                 str(model_tag) + '.h5'
    model_predict_fn_sg = args.model_dir + \
                 str(args.species) + \
                 '_win' + str(args.window_size) + '_' +  \
                 str(args.epochs) +'ep_' + \
                 str(model_tag)
    model_predict_fn = model_predict_fn_sg + '.h5'
    fn_sg = str(args.species) + \
            '_win' + str(args.window_size) + '_' +  \
            str(args.epochs) +'ep_' + \
            str(model_tag)
    model_history_fn = args.output_dir + fn_sg+'_hisroty.csv'
    model_summary_fn = args.output_dir + fn_sg + '_parameters.json'
    return model_predict_fn, model_train_fn, model_history_fn, model_summary_fn


def write_sumery_file(history, summary, history_fn, summary_fn):
    try:
        dict_to_json_write_file(summary, summary_fn)
    except Exception as e:
        raise RuntimeError('！！！[Error] saveing summary hisroty result ', e)
    try:
        pd.DataFrame(history).to_csv(history_fn, index=False)
    except Exception as e:
        raise RuntimeError('！！！[Error] saveing training hisroty result ',e)


def summary_generate(args, start_t):
    end_t = time.time()
    summary={
        'model_type':args.model_type,
        'species':args.species,
        'input_file':args.fast5_base_dir,
        'read_counts':args.read_counts,
        'window_size':args.window_size,
        'epochs':args.epochs,
        'batch_size':args.batch_size,
        'validation_split':args.validation_split,
        'training_time':str(int(end_t-start_t))+' seconds',
    }
    return summary


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
