# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - nanorev_fast5_handeler
 
 @Time: 2019/4/14 4:42 PM
 
 @Author: lotuswang
 
 
"""
import warnings
warnings.filterwarnings('ignore')
import sys
import os
import time
from distutils.version import LooseVersion
import h5py
import numpy as np
import shutil
import logging
from configparser import ConfigParser
from nanorevutils.fileoptions import this_folder
from albacore.config_utils import get_barcoding_options, parse_telemetry_options
from albacore.output_utils import output_summary, write_telemetry
from albacore.input_utils import list_input_files
from albacore import MIN_CFG_VERSION_I, VERSION_I, VERSION
from albacore.pipeline import Pipeline
from albacore.fast5_fastq_data_handler import Fast5FastQDataHandler
from albacore.summary import SummaryCSV
from albacore.log_utils import initialise_logger, log_submission
from albacore.path_utils import initialise_paths
from albacore.config_selector import choose_config
from albacore.telemetry import TelemetryAggregator
from albacore import telemetry_utils
from albacore.aligner_utils import initialise_aligner
from albacore.path_utils import get_default_path


def get_read_data(fast5_fn, basecall_group, basecall_subgroup):
    """
    :param fast5_fn: example:{str} ./input/fast5/id_98490_ch139_read1203_strand.fast5
    :param basecall_group: {str} Default = 'Basecall_1D_000'
    :param basecall_subgroup: {str} Default = 'BaseCalled_template'

    :return:
    **********
    if events.move != 0,save
    if events.move == 2, save 2 elements
    **********
    abs_event_start: event[0]
    start: start times of the bases
    length: bases duration time
    bases:event.model_state[2](move ==1),move == 2
    signal:Raw_signal[abs_event_start:start[-1]+length[-1]]

    all these returns need to be truncated later according to the sam file
    """
    try:
        fast5_data = h5py.File(fast5_fn, 'r')
    except Exception:
        raise NotImplementedError('Error opening file. Likely a corrupted file.')

    try:
        # get albacore version, or if not specified set to 0.0
        albacore_version = LooseVersion(fast5_data['/Analyses/' + basecall_group].attrs['version']
                                        if 'version' in fast5_data['/Analyses/' + basecall_group].attrs
                                        else "0.0")
        if albacore_version <= LooseVersion('0.0'):
            raw_attrs = dict(list(fast5_data['/Raw/Reads/'].values())[0].attrs.items())

            called_dat = fast5_data['/Analyses/' + basecall_group + '/' + basecall_subgroup + '/Events'][()]
            called_dat['start'] = called_dat['start'] * 4000 - raw_attrs['start_time']
            called_dat['length'] = called_dat['length'] * 4000
        else:
            called_dat = fast5_data['/Analyses/' + basecall_group + '/' + basecall_subgroup + '/Events'][()]
    except Exception:
        raise RuntimeError('No events or corrupted events in file. ' +
                           'Likely a segmentation error .')

    start = list()
    bases = list()
    ab_mean = list()
    ab_std = list()
    for p_mean, start_l, p_std, _, model_state, move, _, _ in zip(called_dat['mean'][::-1],
                                                                    called_dat['start'][::-1],
                                                                    called_dat['stdv'][::-1],
                                                                    called_dat['length'][::-1],
                                                                    called_dat['model_state'][::-1],
                                                                    called_dat['move'][::-1],
                                                                    called_dat['p_model_state'][::-1],
                                                                    called_dat['weights'][::-1]):
        model_state = str(model_state)[2:-1]
        start_l = int(start_l)
        if move == 0:
            continue
        elif move == 1:
            start.append(start_l)
            bases.append(model_state[2])
            ab_mean.append(p_mean)
            ab_std.append(p_std)
        elif move == 2:
            start.append(start_l + 2)
            bases.append(model_state[2])
            ab_mean.append(p_mean)
            ab_std.append(p_std)
            start.append(start_l)
            bases.append(model_state[1])
            ab_mean.append(p_mean)
            ab_std.append(p_std)
        else:
            start.append(start_l)
            bases.append(model_state[2])
            ab_mean.append(p_mean)
            ab_std.append(p_std)
    start = start[::-1]
    bases = bases[::-1]
    ab_mean = ab_mean[::-1]
    ab_std = ab_std[::-1]

    try:
        length = np.diff(start)
        length = list(length)
        if start[-1] - start[-2] < 5:
            length.append(3.)
        else:
            length.append(5.)
    except Exception:
        raise RuntimeError('Events is too short or ' +
                           'there are too much zero moves.')

    try:
        read_name = list(fast5_data['/Raw/Reads/'].items())[0][0]
        signal = fast5_data[str('/Raw/Reads/') + str(read_name) + '/Signal/'][()]
    except Exception:
        raise RuntimeError('No signal stored in the file')

    # TODO: channel_info，raw_attrs are necessary or not？
    # fast5_info = fast5_data['UniqueGlobalKey/channel_id'].attrs
    # raw_attrs = dict(fast5_data['/Raw/Reads/'].values()[0].attrs.items())

    fast5_data.close()
    if len(signal) < int(start[-1] + length[-1]):
        raise RuntimeError('Signal is shorter than the Events')
    else:
        abs_event_start = start[0]
        start = np.array(start) - abs_event_start
        length = np.array(length)
        # if len(quality)!=len(start):
        #     print(fast5_fn,'quality error')
        return abs_event_start, start, length, bases, signal, ab_mean, ab_std


class Opt(object):
    def __init__(self, default_path, fast5_base_dir, temp_dir):
        self.default_path = default_path
        self.input = fast5_base_dir
        self.save_path = temp_dir
        self.worker_threads = 1
        self.resume = ''
        self.flowcell = 'FLO-MIN106'
        self.kit = 'SQK-LSK108'
        self.debug = False
        self.recursive = False
        self.files_per_batch_folder = 0
        self.output_format = 'fastq'
        self.reads_per_fastq_batch = 1
        self.disable_filtering = True
        self.disable_pings = True
        self.version = '1.0'
        self.logfile = 'pipeline.log'
        self.summfile = 'sequencing_summary.txt'
        self.cfgfile = 'configuration.cfg'
        self.data_path = None
        self.barcoding = False
        self.align_ref = None
        self.config = None


def handle_opts(default_path, fast5_base_dir, temp_dir):
    result_opts = Opt(default_path, fast5_base_dir, temp_dir)
    return result_opts


def opt_main(opts):

    opts.read_path = initialise_paths(opts.save_path)
    logger = initialise_logger(opts.save_path, opts.debug, logfilename=opts.logfile, resume_mode=opts.resume)
    if opts.data_path is None:
        opts.data_path = opts.default_path

    if opts.config is None:
        opts.config, barcoding = choose_config(opts.data_path, opts.flowcell, opts.kit)
        if barcoding: # If the flowcell+kit includes barcoding we always run barcoding
            opts.barcoding = True
            logger.info('Kit {} includes barcodes "--barcoding" option set to True'.format(opts.kit))
    else: # Otherwise we look for it first as an abs path, then in the data path
        if not os.path.exists(opts.config):
            opts.config = os.path.join(opts.data_path, opts.config)
            if not os.path.exists(opts.config):
                raise Exception('Config file "{}" does not exist.'.format(opts.config))

    config = ConfigParser(interpolation=None)
    config.read(opts.config)

    # Clear out any compatibility information in the config file
    config.remove_section('compatibility')

    # The python loader expects a model path to be set, but we don't need it elsewhere.
    if config.has_section('basecaller') and not config.has_option('basecaller', 'model_path'):
        config.set('basecaller', 'model_path', opts.data_path)
    elif config.has_section('basecaller_hmm') and not config.has_option('basecaller_hmm', 'model_path'):
        config.set('basecaller_hmm', 'model_path', opts.data_path)
    elif config.has_section('basecaller_ocl') and not config.has_option('basecaller_ocl', 'model_path'):
        config.set('basecaller_ocl', 'model_path', opts.data_path)

    # in case a specific version is specified in the config file, check if this
    # version is in the range [min_config_version, current_version]
    config_version = config.get('pipeline', 'albacore_version', fallback=None)
    if config_version:
        config_version_info = tuple([int(num) for num in
                                     config_version.split('.')])
        if (config_version_info < MIN_CFG_VERSION_I
            or config_version_info > VERSION_I):
            msg = ("Current version of albacore ({}) is incompatible with the "
                   "config file's version ({}).".format(VERSION,
                                                        config_version))
            logger.error(msg)
            raise Exception(msg)
        config.remove_option('pipeline', 'albacore_version')

    get_barcoding_options(config, opts.barcoding, opts.data_path, logger,
                          opts.resume)
    # Now the config file is complete, all "usual" command-line arguments have
    # been parsed, and all sub-config files have been added. Now look for command-line overrides:

    # find external aligners and references; delete respective config section
    # if alignment or calibration strand detection is disabled
    references_loaded = \
        initialise_aligner(config, opts.align_ref, opts.data_path, 'aligner') + \
        initialise_aligner(config, None, opts.data_path, 'calib_detector')

    basecall_type = config.get('pipeline', 'basecall_type')
    raw_basecalling = not config.has_section('event_detector')
    ocl_basecall = config.has_section('basecaller_ocl')
    basecaller_section = 'basecaller_ocl' if ocl_basecall else 'basecaller'
    if basecall_type != 'linear' and raw_basecalling:
        raise Exception('Basecalling without event detection (a.k.a. raw basecalling) is currently '
                        'only available for 1D.')
    desc_base = {'linear': '1d',
                 'full_2d': '2d'}
    desc_file = 'layout_' + ('raw_' if raw_basecalling else '') + 'basecall_' \
        + desc_base[basecall_type] + ('_ocl' if ocl_basecall else '') + '.jsn'
    logger.info("Using config file: {}".format(desc_file))
    config.set('pipeline', 'desc_file', os.path.join(opts.data_path, desc_file))

    if config.has_option(basecaller_section, 'opencl_kernel_path'):
        config.set(basecaller_section, 'opencl_kernel_path',
                   os.path.join(opts.data_path,
                                config.get(basecaller_section, 'opencl_kernel_path')))

    files, additional_summary_lines = list_input_files(opts) # will set opts.resume to False
                                                             # if no summary file is present

    # Initialize the summary file handler.
    include_calib = bool(config.get('calib_detector', 'method'))
    include_alignment = bool(config.get('aligner', 'method'))
    include_barcoding = config.has_option('barcode_detector',
                                          'arrangements_files')
    logger.info('include calibration strand detection: '
                '{}'.format(include_calib))
    logger.info('include alignment: {}'.format(include_alignment))
    logger.info('include barcoding: {}'.format(include_barcoding))

    telemetry_enabled = not opts.disable_pings
    (telemetry_urls,
     telemetry_segment_duration,
     telemetry_analysis_id) = parse_telemetry_options(config)
    telemetry_output_file = os.path.join(opts.save_path,
                                         'sequencing_telemetry.js')
    telemetry_aggregator = TelemetryAggregator(
        telemetry_segment_duration,
        albacore_opts=opts,
        include_1d_basecalling=True,
        include_1dsq_basecalling=False,
        include_calibration=include_calib,
        include_alignment=include_alignment,
        include_barcoding=include_barcoding,
        analysis_id=telemetry_analysis_id)
    logger.info('telemetry enabled: {}'.format(telemetry_enabled))
    summary_fields = output_summary(basecall_type, include_calib,
                                    include_barcoding, include_alignment)
    try:
        with SummaryCSV(os.path.join(opts.save_path, opts.summfile),
                        summary_fields,
                        start_new_file=not opts.resume) as summary:
            for filename in additional_summary_lines:
                summary.write({'filename': os.path.basename(filename)})
            process_pipeline(opts, config, files, summary, telemetry_aggregator)
        logger.info('Done.')
    except Exception as exc:
        logger.exception('Fatal error in pipeline')
        print('Fatal error in pipeline:', exc, file=sys.stderr)
        raise
    logger.info('Writing telemetry')
    write_telemetry(telemetry_aggregator, telemetry_output_file)
    if telemetry_enabled:
        for url in telemetry_urls:
            logger.info('Pinging telemetry to {}'.format(url))
            telemetry_utils.ping(telemetry_aggregator, url)
    if references_loaded:
        unload_commands = [config.get('aligner', 'unload_command').
                           format(executable=config.get('aligner',
                                                        'executable'),
                                  reference=ref)
                           for ref in references_loaded]
        print('\n    CAVE: Reference is still loaded. To unload, call:')
        for cmd in set(unload_commands):
            print('    ' + cmd)
        print('    ')
    return 0


def process_pipeline(opts, config, files, summary, telemetry_aggregator=None):
    logger = logging.getLogger('albacore')

    # Zero threads actually means turning off multithreading. But for purposes
    # of determining whether the pipeline is finished working, we need to know
    # the actual number of workers, which is still one.
    num_workers = opts.worker_threads if opts.worker_threads > 0 else 1
    # channels = list(range(1, num_workers + 1))
    reverse_direction = config.getboolean('basecaller', 'reverse_direction',
                                          fallback=False)
    u_substitution = config.getboolean('basecaller', 'u_substitution',
                                       fallback=False)
    new_config = os.path.join(opts.save_path, opts.cfgfile)
    with open(new_config, 'w') as cfg_out:
        config.write(cfg_out)

    # Initialize the basecall pipeline.
    data_handler = Fast5FastQDataHandler(new_config, opts)
    # Note that we pass opts.worker_threads (which may be zero), not num_workers.
    pipeline = Pipeline(opts.worker_threads, data_handler, new_config,
                        reverse_direction=reverse_direction,
                        u_substitution=u_substitution,
                        telemetry_aggregator=telemetry_aggregator)

    num_files = len(files)
    file_index = 0
    num_complete = 0

    # only copy files if Fast5 output is enabled
    if 'fast5' in data_handler.output_format:
        copied_files = []
        num_copied = 0
    else:
        copied_files = files
        num_copied = num_files
    no_more_files = False

    # pbar = get_progress_bar(opts.debug, num_files)

    PASSING_DATA = 0
    FINISHING_JOBS = 1
    ALL_JOBS_COMPLETE = 2

    pipeline_state = PASSING_DATA
    while pipeline_state != ALL_JOBS_COMPLETE:
        if pipeline_state == FINISHING_JOBS:
            num_cached_reads = pipeline.finish_all_jobs()
            pipeline_state = ALL_JOBS_COMPLETE
        else:
            num_cached_reads = pipeline.cache_completed_calls()
            free_workers = pipeline.workers_ready()
            if free_workers == num_workers and no_more_files:
                # All workers are idle, and there is no more data to process. We're done.
                pipeline_state = FINISHING_JOBS
            will_sleep = True
        if free_workers == 0 and num_cached_reads == 0:
            if num_copied < num_files:
                try:
                    fname = files[num_copied]
                    destination = os.path.join(opts.read_path,
                                               os.path.basename(fname))
                    subpath = os.path.dirname(destination)
                    if not os.path.exists(subpath):
                        os.makedirs(subpath)
                    shutil.copyfile(fname, destination)
                    copied_files.append(destination)
                    num_copied += 1
                    will_sleep = False
                except Exception:
                    logger.error('Could not access file {}. '
                                 'File will be skipped.'.format(fname))
                    files.remove(fname)
                    num_files = len(files)

        if free_workers > 0 and file_index < num_files:
            # Submit new jobs before processing complete ones, so that our worker
            # threads are not left idle.
            submit = min(free_workers, num_files - file_index)
            for i in range(submit):
                if num_copied <= file_index:
                    try:
                        fname = files[file_index]
                        destination = os.path.join(opts.read_path,
                                                   os.path.basename(fname))
                        subpath = os.path.dirname(destination)
                        if not os.path.exists(subpath):
                            os.makedirs(subpath)
                        shutil.copyfile(fname, destination)
                        copied_files.append(destination)
                        num_copied += 1
                    except Exception:
                        logger.error('Could not access file {}. '
                                     'File will be skipped.'.format(fname))
                        files.remove(fname)
                        num_files = len(files)

                if num_copied > file_index:
                    short_name = os.path.basename(copied_files[file_index])
                    logger.info(log_submission(short_name))
                    try:
                        status = pipeline.submit_read(copied_files[file_index])
                        if not status:
                            logger.warning('Could not extract read data from '
                                           'file "{}".'.format(short_name))
                    except Exception as e:
                        logger.exception('Could not submit file "{}, '
                                         'error: {}".'.format(short_name, e))
                    file_index += 1
                    will_sleep = False

        if file_index == num_files:
            no_more_files = True

        if num_cached_reads > 0:
            # Process the data from the completed workers.
            summary_data = pipeline.process_cached_calls()
            for line in summary_data:
                line['filename'] = '{}.fast5'.format(line['label'])
                num_complete += 1
                if 'error_message' in line:
                    logger.warning('Error processing '
                                   'file "{}".'.format(line['filename']))
                    logger.warning(line['error_message'])
                else:
                    summary.write(line)
                    logger.info('Finished processing '
                                'file "{}".'.format(line['filename']))
            will_sleep = False

        if will_sleep:
            # No workers are ready, and we haven't done any work. Take a nap.
            time.sleep(0.01)