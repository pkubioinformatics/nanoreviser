
# NanoReviser: An Error-correction Tool for Nanopore Sequencing Based on a Deep Learning Algorithm

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Version](#version)
* [Installation](#installation)
* [Usage](#usage)
    * [NanoReviser.py](#NanoReviser.py)
    * [NanoReviser_train.py](#NanoReviser_train.py)
* [Examples](#examples)
    * [For basecalling revising](#for-basecalling-revising)
    * [For training NanoReviser](#for-training-nanoreviser)
* [Citation](#citation)
* [Contact](#contact)

## Introduction

Nanopore sequencing is regarded as one of the most promising third-generation sequencing (TGS) technologies. However, the nanopore sequencing reads are susceptible to a fairly high error rate owing to the difficulty in identifying the DNA bases from the complex electrical signals. Here we proposed a DNA basecalling reviser, **NanoReviser**, based on a deep learning algorithm to correct the basecalling errors introduced by current basecallers provided by default. In our module, we re-segmented the raw electrical signals based on the basecalled sequences provided by the default basecallers. By employing convolution neural networks (CNNs) and bidirectional long short-term memory (Bi-LSTM) networks, we took advantage of the information from the raw electrical signals and the basecalled sequences from the basecallers. 

**NanoReviser**, as a post-basecalling reviser, significantly improving the basecalling quality. The NanoReviser package is freely available at https://github.com/pkubioinformatics/NanoReviser.


## Requirements

+ [Anaconda](https://www.anaconda.com/)
+ [h5py, version >= 1.8.12](http://www.hdfgroup.org/HDF5/)
+ [Python version >= 3.6.9](https://www.python.org/)
+ [numpy version >= 1.17.3](http://www.numpy.org/)
+ [pandas, version >= 0.25.0](http://pandas.pydata.org/)
+ [TenserFlow, version >= 1.12.0](https://www.tensorflow.org/)
+ [Keras, version 1.2.2](https://https://github.com/keras-team/keras/)
+ [GraphMap, version >= 0.5.2](https://github.com/isovic/graphmap/)

## Version

+ NanoReviser 1.0 (Tested on MacOS 10.11 and Linux_64, including CentOS 7 and Ubuntu 18.04)


## Installation


#### Install NanoReviser using git

clone NanoReviser package

    $ git clone https://github.com/pkubioinformatics/nanoreviser.git
    
change directory to NanoReviser

    $ cd nanoreviser

If you currently have TensorFlow installed on your system, we would advise you to create a virtual environment to install NanoReviser. If you want to do so, we recommend the user-friendly [Anaconda](https://www.anaconda.com/).

You will create an eviroment named nanorev for NanoReviser and install all dependencies through conda

**Note : you need to replace the path where you installed anaconda**

For Linux Ubuntu 18.04 system:

     $sudo apt-get install libhdf5-serial-dev hdf5-tools

For Linix Centos 7 system:
    
     $ sudo yum install https://download-ib01.fedoraproject.org/pub/epel/7/x86_64/Packages/h/hdf5-1.8.12-11.el7.x86_64.rpm


For linux with gpu 

      $ conda env create -n nanorev /**Your_Path_to_Anaconda**/envs/nanorev/ -f ./enviroment/NanoReviser.yaml 
	  $ conda activate nanorev
 

For linux just with cpu
	
	  $ conda env create -n nanorev /**Your_Path_to_Anaconda**/envs/nanorev/ -f ./enviroment/NanoReviser_cpu.yaml 
	  $ conda activate nanorev
	  $ conda install tensorflow==1.12.0

For macOS with cpu 

      $ conda env create -n nanorev /**Your_Path_to_Anaconda**/envs/nanorev/ -f ./enviroment/NanoReviser_macOS.yaml 
      $ conda activate nanorev

Please run the unitest to make sure NanoReviser installed properly.

      $ conda activate nanorev
      $ sh unitest.sh

If both NanoReviser.py and NanoReviser_train.py are available, you will get

       Congratulations, please have fun with NanoReviser :)
    

## Usage


#### NanoReviser.py

An ONT basecalling reviser based on deep learning

    usage:
           python NanoReviser.py [options] -d <fast5_files> -o <output_path> -f <output_format, default=fasta>

           [example]
           python NanoReviser.py -d ./unitest/test_data/fast5/ -o ./unitest/nanorev_output/ -F fasta -S ecoli

	  usage: 
           python NanoReviser.py [options]

	  An ONT basecalling reviser based on deep learning
	
    Options:
    --version                                              show program's version number and exit
    -h, --help                                             show this help message and exit
    -d FAST5_BASE_DIR, --fast5_base_dir=FAST5_BASE_DIR     path to the fast5 files
    -o OUTPUT_DIR, --output_dir=OUTPUT_DIR                 path to store the output files
    -S SPECIES, --species=SPECIES                          species model to load which located in ./model/, 
                                                           default is "human"
    -F OUTPUT_FORMAT, --output_format=OUTPUT_FORMAT        format of the output files, default is fasta
    --thread=THREAD                                        thread, default is 100
    -t TEMP_DIR, --tmp_dir=TEMP_DIR                        path to the tmp dir, which is used to store the 
                                                           preprocessing files
    -e FAILED_READS_FILENAME                               document to log the failed reads, default is
                                                           failed_read.txt
    


#### NanoReviser_train.py

A training tools for generation model files for NanoReviser

    usage:
           python NanoReviser_train.py [options] 

	  An ONT basecalling reviser based on deep learning

	  Options:
    --version                                              show program's version number and exit
    -h, --help                                             show this help message and exit
    -d FAST5_BASE_DIR, --fast5_base_dir=FAST5_BASE_DIR     path to the fast5 files
    -r REFERENCE, --reference=REFERENCE                    reference genome for labeling the training data
    -m OUTPUT_MODEL, --output_dir=OUTPUT_MODEL             name of the dir to store model1 and model2
    -o OUTPUT_DIR, --output_dir=OUTPUT_DIR                 path to store the output summery files

    -b BATCH_SIZE, --batch_size=BATCH_SIZE                 batch size of trainig NanoReviser, default is 256
    -e EPOCHS, --epochs=EPOCHS                             epochs of training NanoReviser, defualt is 50
    -w WINDOW_SIZE, --window_size=WINDOW_SIZE              window size for slicing the read input, defualt is 13
    -c READ_COUNT, --read_count=READ_COUNT                 the number of read included in the training data, must 
                                                           smaller than the number of files stored in fast5_base_dir, 0 for use all the files in the fast5_base_dir and defult is 0.
    --validation_split                                     validation data size, default is 0.01, which means 1% 
                                                           reads in fast5_base_dir would be used as validation data.

    --thread=THREAD                                        thread, default is 1
    --model_type=MODEL_TYPE                                'both', 'model1' or 'model2', default is 'both'
    --mapper_exe=MAPPER_EXE                                the align tool for generate the lable of training 
                                                           data, default is 'graphmap'


## Example

#### For basecalling revising

**NanoReviser.py** : An ONT basecalling reviser based on deep learning

For revising the fast5 files in ./unitest/test_data/fast5/ in order to get fasta files,the command line would be:

    $ conda activate nanorev  #activate the python enviroment for nanoreviser
    $ pyton NanoReviser.py -d ./unitest/test_data/fast5/ -o ./unitest_nanorev_results/ -S ecoli -F fasta

For revising the fast5 files in ./unitest/test_data/fast5/ in order to get fastq files,the command line would be:

    $ conda activate nanorev  #activate the python enviroment for nanoreviser
    $ pyton NanoReviser.py -d ./unitest/test_data/fast5/ -o ./unitest_nanorev_results/ -S ecoli -F fastq

Please run the following command in oder to get the entire fasta or fastq file contains all reads in fasta5's dir:

    $ cat ./nanorev_results/*.fasta > nanorev_results.fasta 
or
    
    $ cat ./nanorev_results/*.fastq > nanorev_results.fastq



#### For training NanoReviser

**NanoReviser_train.py** A training tools for generation model files for NanoReviser

For training NanoReviser by data in ./unitest/training_data/fast5/ and reference genome in ./unitest/training_data/reference.fasta in order to get model files in ./model/unitest/ and result files in ./unitest/training_result/,the command line would be:

    $ conda activate nanorev  #activate the python enviroment for nanoreviser
    $ pyton NanoReviser_train.py -d ./unitest/training_data/fast5/ -r ./unitest/training_data/reference.fasta -o ./unitest_training_results/ -S unitest

This command will generate two model files in ./model/unitest and for summery files in ./unitest_training_data/.

Please note that the training process of NanoReviser_train.py could take quite a long time, we highly recommend to use screen command to run NanoReviser_train.py as follow:

    $ screen -S nanorev_train
    $ conda activate nanorev  #activate the python enviroment for nanoreviser
    $ pyton NanoReviser_train.py -d ./unitest/training_data/fast5/ -r ./unitest/training_data/reference.fasta -o ./unitest_training_results/ -S unitest -b 256 -w 13 -e 50 -c 0 --validation_slipt=0.01 --model_type=both

There would be four result files in ./unitest_training_results/:

  (1)  **unitest_win13_50ep_model1_history.csv** (Training history of NanoReviser model1)

  (2)  **unitest_win13_50ep_model1_summery.json** (Paramaters of of NanoReviser model1)

  (3)  **unitest_win13_50ep_model2_history.csv** (Training history of NanoReviser model2)

  (4)  **unitest_win13_50ep_model2_summery.json** (Paramaters of of NanoReviser model2)

And two model files in ./model/unitest:

  (1)  **unitest_win13_50ep_model1.h5** 

  (2)  **unitest_win13_50ep_model2.h5** 


## Citation


Luotong Wang, Li Qu, Longshu Yang, Yiying Wang and Huaiqiu Zhu\*; NanoReviser: An Error-correction Tool for Nanopore Sequencing Based on a Deep Learning Algorithm


## Contact


Please direct your questions to: Dr. Huaiqiu Zhu, [hqzhu@pku.edu.cn](hqzhu@pku.edu.cn)


