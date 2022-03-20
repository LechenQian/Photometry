# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 17:24:30 2022

@author: qianl
"""
#%% load function and package


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
from sklearn.linear_model import LinearRegression

#%% load data
mouse_name = 'FgDA_01'
gcamp_data = load_pickleddata('D:/PhD/Photometry/DATA/round_20220307/processed/corrected_gcamp_data.pickle')
mouse_data = load_pickleddata('D:/PhD/Photometry/DATA/round_20220307/processed/{}_stats.pickle'.format(mouse_name))
#%%
session= 2
if mouse_name == 'FgDA_01':
    
    good sites = ['AMOT','PMOT','ALOT','PLOT','MNacS','LNacS']
else:
    good_sites = mouse_data['df_bpod_doric'][str(session)+'_'+mouse_data['all_days'][session]]

num_neuron = len(good_sites)
rel_contrib_mat = np.zeros([num_neuron,5]) # create a empty matrix to store the relative contribution of each variable

GLM_DICT = {}
for index,site in enumerate(good_sites):
    beh_df = mouse_data['df_bpod_doric'][str(session)+'_'+mouse_data['all_days'][session]]['dataframe']
    
    binned_licking = []
    binned_dff = []
    diff_binned_licking = []
    trialnum = []
    water = np.zeros([7,7]) ？？？
    go_odor = np.zeros([7,7])
    nogo_odor = np.zeros([7,7])
    
    GLM_df = pd.DataFrame()
    # fluo_signal = gcamp_data[mouse_name][TT]['signal_gcamp']
    for i, row in beh_df.iterrows():
        count_go = 0
        count_nogo = 0
        count_omit = 0
        if row['Trialtype'] != 'background':
            
            
            binned_licking += np.histogram(row['lickings'][0],bins = 45,range = (0,9))[0]/(1/5)

            diff_licking = np.diff(bin_licking(row_df=row,trialtype=row['trialtype'],binsize = 1))
            diff_licking = np.insert(diff_licking,0,0)
            diff_binned_licking += list(diff_licking)
            
            # whole-trial variable
            trialnum += [i]*45
            # events variable
            x = np.linspace(0,9,45)
            knots = sb.get_knots_bsplines(x,df = 7,)#spacing = 'equal')
    
            basis = sb._eval_bspline_basis(x,degree = 3,knots=knots)[0]
            temp_water = np.zeros([45,7])
            temp_go_odor = np.zeros([45,7])
            temp_nogo_odor = np.zeros([45,7])
            if row['trialtype'] == 'go':
                
                binned_dff += gcamp_data[mouse_name][TT]['signal_gcamp'][site][count_go,:,session] # replace nan value with zero and convert to 45 element
                water_event = np.zeros([1,45])
                water_event[27] = 1
                go_odor_event = np.zeros([1,45])
                go_odor_event[9:14] = np.ones([1,5])
                ## 今天写到这里写不动了
                for i in range(7):
                    conv_water = np.convolve(basis[:,i],water_event,mode = 'full')
                    temp_water[i,:] = conv_water[0:7]
                    conv_go = np.convolve(basis[:,i],go_odor_event,mode = 'full')
                    temp_go_odor[i,:] = conv_go[0:7]
                count_go += 1
            elif row['trialtype'] == 'go_omit':
                binned_dff += temp
                go_odor_event = np.array([0,1,0,0,0,0,0])
                for i in range(7):
                    conv_go = np.convolve(basis[:,i],go_odor_event,mode = 'full')
                    temp_go_odor[i,:] = conv_go[0:7]
                count_nogo += 1
            elif row['trialtype'] == 'no_go':
                binned_dff += temp
                water_event = np.array([0,0,0,0,0,0,0])
                nogo_odor_event = np.array([0,1,0,0,0,0,0])
                for i in range(7):
                    conv_water = np.convolve(basis[:,i],water_event,mode = 'full')
                    temp_water[i,:] = conv_water[0:7]
                    conv_nogo = np.convolve(basis[:,i],nogo_odor_event,mode = 'full')
                    temp_nogo_odor[i,:] = conv_nogo[0:7]
                count_omit += 1
            water = np.hstack((water,temp_water))
            go_odor = np.hstack((go_odor,temp_go_odor))
            nogo_odor = np.hstack((nogo_odor,temp_nogo_odor))
    
    water = water[:,7:]
    go_odor = go_odor[:,7:]
    nogo_odor = nogo_odor[:,7:]
    GLM_df = pd.DataFrame({'NeuroSignal_dff':binned_dff,'licking':binned_licking,'trialnum':trialnum,'diff_licking':diff_binned_licking,'square_licking':np.power(binned_licking,2)})
    
    for i in range(7):
        GLM_df['water{}'.format(i)] = water[i,:]
    for i in range(7):
        GLM_df['go_odor{}'.format(i)] = go_odor[i,:]    
    for i in range(7):
        GLM_df['nogo_odor{}'.format(i)] = nogo_odor[i,:]    
    GLM_df.describe()