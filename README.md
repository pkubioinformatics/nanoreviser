
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
+ [Python version >= 3.6.9](https://www.python.org/)
+ [numpy version >= 1.17.3](http://www.numpy.org/)
+ [pandas, version >= 0.25.0](http://pandas.pydata.org/)
+ [TenserFlow, version >= 1.12.0](https://www.tensorflow.org/)
+ [Keras, version 1.2.2](https://https://github.com/keras-team/keras/)
+ [h5py, version >= 2.5.0](http://www.h5py.org/)
+ [GraphMap, version >= 0.5.2](https://github.com/isovic/graphmap/)

## Version

+ NanoReviser 1.0 (Tested on Linux_64, including CentOS 6.5 and Ubuntu 16.04)


## Installation


#### Install NanoReviser using git

clone NanoReviser package

    $ git clone https://github.com/pkubioinformatics/nanoreviser.git
    
change directory to NanoReviser

    $ cd nanoreviser

If you currently have TensorFlow installed on your system, we would advise you to create a virtual environment to install NanoReviser. If you want to do so, we recommend the user-friendly [Anaconda](https://www.anaconda.com/).

You will create an eviroment named nanorev for NanoReviser and install all dependencies through conda
 
 1. for linux with gpu 

	#you need to replace the path where you installed anaconda, the default path would be ~/anaconda3/

	$ conda env create -n nanorev /**Your_Path_to_Anaconda**/envs/nanorev/ -f NanoReviser.yaml  
	$ conda activate nanorev
 
 2. for linux just with cpu

	#you need to replace the path where you installed anaconda, the default path would be ~/anaconda3/

	$ conda env create -n nanorev_cpu /**Your_Path_to_Anaconda**/envs/nanorev/ -f NanoReviser_cpu.yaml  
	$ conda activate nanorev_cpu	
 
 3. for macOS
    #you need to replace the path where you installed anaconda, the default path would be ~/anaconda3/
    
    $ conda env create -n nanorev /**Your_Path_to_Anaconda**/envs/nanorev/ -f NanoReviser_macOS.yaml  
	$ conda activate nanorev

Please run the unitest in oder to make sure NanoReviser installed properly.

    $ sh unitest.sh
    

## Usage


#### NanoReviser.py

An ONT basecalling reviser based on deep learning

    usage:
           python NanoReviser.py [options] -d <fast5_files> -o <output_path>

	usage: 
           python NanoReviser.py [options]

	An ONT basecalling reviser based on deep learning

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit
	  -d, --fast5_base_dir
	  -o, --output_dir

#### NanoReviser_train.py

A training tools for generation model files for NanoReviser

    usage:
           python NanoReviser_train.py [options] -d <fast5_files> -m <output_model>

	usage: 
           python NanoReviser_train.py [options]

	An ONT basecalling reviser based on deep learning

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit
	  -d, --fast5_base_dir
	  -o, --output_dir

### Authour

Luotong Wang, Li Qu, Longshu Yang, Yiying Wang Huaiqiu Zhu; NanoReviser: An Error-correction Tool for Nanopore Sequencing Based on a Deep Learning Algorithm


### Contact

Please direct your questions to: Dr. Huaiqiu Zhu, [hqzhu@pku.edu.cn](hqzhu@pku.edu.cn)


