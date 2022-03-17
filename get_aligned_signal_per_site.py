# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 18:06:58 2022

@author: qianl
"""

from datetime import datetime
from datetime import date
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
import os
import random
import matplotlib as mpl
import re
import csv
import pickle
import sys
# sys.path.insert(0,os.path.join(path,'functions'))
from load_all import Mouse_data
from load_all import pickle_dict
from load_all import load_pickleddata


mouse_names = ['FgDA_01','FrgD1_01','FrgD2_01','FrgD2_02']
path = 'D:/PhD/Photometry/DATA/round_20220307'
for mouse_id in mouse_names:    
    
    load_path = os.path.join(path,'processed/{0}_stats.pickle'.format(mouse_id))
    mouse = load_pickleddata(load_path)
    
    #event plot with trials and iscorerct data
    
    # assign two df 
    mouse_trials = mouse.df_bpod_doric
    
    # choose a date
    all_days = mouse.all_days