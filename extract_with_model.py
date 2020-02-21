# -*- coding: utf-8 -*-
"""
 @File: correcter - extract_direct
 
 @Time: 2020/2/8 2:10 PM
 
 @Author: lotuswang
 
 
"""
import os
import time
from optparse import OptionParser
from multiprocessing import Pool
from deepreutils import *
from pathlib import Path


os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

def get_args():
    optParser = OptionParser(usage="%prog [-d] [-r] [-o]", version="%prog 1.0")

    optParser.add_option('-d', '--fast5_base_dir', action='store', type="string",
                         dest='fast5_base_dir', default='../MNT_test_set/fast5/')
    optParser.add_option('-t', '--tmp_dir', action='store', type="string", dest='temp_dir',
                         default='./tmp/')
    optParser.add_option('-o', '--output_dir', action='store', type="string", dest='output_dir',
                         default='./nanorev_output/')
    optParser.add_option('-r', '--reference_genome', action='store', type="string", dest='genome_fn',
                         default='../MNT_test_set/NC_000913.3_ecoli_k12_MG1655.fasta')
    optParser.add_option("--thread", action="store", type="int", dest="thread",
                         default=100)
    optParser.add_option('-m', '--mapper_dir', action='store', type="string", dest='graphmap_exe',
                         default='graphmap')
    optParser.add_option('-e', '--failed_read', action='store', type="string", dest='failed_reads_filename',
                         default='failed_reads.txt')
    optParser.add_option('-g', '--basecall_group', action='store', type="string", dest='basecall_group',
                         default='Basecall_1D_000')
    optParser.add_option('-s', '--basecall_subgroup', action='store', type="string", dest='basecall_subgroup',
                         default='BaseCalled_template')
    optParser.add_option('-L', '--output_format', action='store', type="string", dest='output_format',
                         default='sam')
    optParser.add_option('--model1_predict_dir', action='store', type="string", dest='model1_predict_dir',
                         default='../model/MNT_low_local_att_model1_win11_predict_50_10ep.h5')
    optParser.add_option('--model2_predict_dir', action='store', type="string", dest='model2_predict_dir',
                         default='../model/MNT_low_local_att_model2_win11_predict_50_10ep.h5')
    optParser.add_option("-v", "--virsion", action="store_true", dest="virsion",
                         help="version of NanoReviser")




    (tmp_args, _) = optParser.parse_args()
    if tmp_args.virsion:
        print("The virsion of NanoReviser : 0.1 ")
    return tmp_args

def check_path(my_path):
    my_file = Path(str(my_path))
    if not my_file.exists():
        try:
            os.mkdir(my_path)
        except Exception as e:
            print('！！！[Error] make dir ' + my_file + e)


# 存在
def provide_fasta(name, fast5_fn_sg, args):

    print('Run task %s (%s)...' % (name, os.getpid()))
    fast5_fn = os.path.join(args.fast5_base_dir, fast5_fn_sg)
    try:
        (read_start_rel_to_raw, starts_rel_to_read, event_length, event_bases,
         or_raw_signal, raw_mean, raw_std) = get_read_data(fast5_fn, args.basecall_group, args.basecall_subgroup)
    except Exception as e:
        print('！！！[Error] fast5 file: ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        signal = or_raw_signal[int(read_start_rel_to_raw):]
        signal_list, signal_mean, signal_std, shift, scale = signal_segmentation(signal,
                                                                                 starts_rel_to_read,
                                                                                 int(event_length[-1]))
        signal_mean = signal_mean / shift
        signal_std = signal_std / scale
    except Exception as e:
        print('[！！！Error] signal: ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        read_bases = np.array(pd.Series(event_bases).apply(get_base_color)) / 300.0
        clean_event_length = np.array(event_length) / 10.0
        x = np.vstack([read_bases, signal_mean, signal_std, clean_event_length, raw_mean, raw_std])
        x = x.T
        x = np.array(x, dtype='float')
        x_tmp = list()
        for i in range(len(x) - 11):
            tmp_x = x[i:i + 11]
            x_tmp.append(tmp_x)
        x = np.array(x_tmp)
        del x_tmp
        signal_x_tmp = list()
        for i in range(len(signal_list) - 11):
            tmp_x = signal_list[i:i + 11]
            signal_x_tmp.append(tmp_x)
        # x_train = np.array(tmp_x)p
        print(len(signal_x_tmp))
        signal_x = np.array(signal_x_tmp)
    except Exception as e:
        print('[！！！Error] input features: ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        att_model1_predict = get_model1()
        att_model2_predict = get_model2()
        att_model1_predict.load_weights(args.model1_predict_dir)
        att_model2_predict.load_weights(args.model2_predict_dir)
        y_pred_total1 = att_model1_predict.predict([signal_x[:, :, :, np.newaxis], x])
        y_pred_total2 = att_model2_predict.predict([signal_x[:, :, :, np.newaxis], x])
        y_pred = np.argmax(y_pred_total1, axis=1)
        y_pred2 = np.argmax(y_pred_total2, axis=1)
    except Exception as e:
        print('[！！！Error] prediction : ' + fast5_fn_sg.split('.')[0] + str(e))
    try:
        y_read = get_base_1(event_bases, y_pred, y_pred2)
        out_fasta_fn = args.output_dir + fast5_fn_sg.split('.')[0] + '_out.fasta'
        if not os.path.exists(args.output_dir):
            os.mkdir(args.output_dir)
        prep_read_fasta(fast5_fn, out_fasta_fn, list(y_read))
        print('[p:::] ' + fast5_fn_sg.split('.')[0] + '_out.fasta was saved......')
    except Exception as e:
        print('[！！！Error] stroring : ' + fast5_fn_sg.split('.')[0] + str(e))

if __name__ == '__main__':
    ar_args=get_args()
    print(ar_args)
    # genome_fn = ar_args.genome_fn
    # genome_index = parse_fasta(genome_fn)
    # print(genome_fn, 'has been load......')
    check_path(ar_args.temp_dir)
    check_path(ar_args.output_dir)

    fast5_fns = os.listdir(ar_args.fast5_base_dir)[:2]
    if int(ar_args.thread) > len(fast5_fns):
        pool_size = len(fast5_fns)
    else:
        pool_size = int(ar_args.thread)


    # 存在
    # run_time = int(len(fast5_fns)/pool_size)
    run_time = int(len(fast5_fns)/pool_size)
    start_time = time.time()
    p = Pool(pool_size)
    for j in range(run_time):
        for i in range(pool_size):
            p.apply_async(provide_fasta,  args=(i, fast5_fns[pool_size*j+i],ar_args,))
    print('Waiting for all subprocesses done...')
    p.close()
    p.join()
    print('All subprocesses done.')
    # for fast5_fn_sg in fast5_fns:
    #     provide_fasta(fast5_fn_sg, fast5_fn_sg, ar_args,genome_index)
    end_time = time.time()
    print(print("NanoReviser time consuming:%.2f seconds" % (end_time - start_time)))

    try:
        os.removedirs(ar_args.output_dir)
    except Exception as e:
        print(e)



