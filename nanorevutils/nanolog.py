# -*- coding: utf-8 -*-
"""
 @File: nanoreviser - nanolog
 
 @Time: 2020/4/15 9:27 PM
 
 @Author: lotuswang
 
 
"""

import logging

def logger_config(log_path,logging_name):
    logger = logging.getLogger(logging_name)
    logger.setLevel(level=logging.DEBUG)
    handler = logging.FileHandler(log_path, encoding='UTF-8')
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    console = logging.StreamHandler()
    console.setLevel(logging.DEBUG)
    logger.addHandler(handler)
    logger.addHandler(console)
    return logger