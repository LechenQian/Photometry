# -*- coding: utf-8 -*-
"""
Created on Mon Oct  3 21:21:17 2022

@author: qianl
analysis between regions
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
from b_load_processed_mat_save_pickle_per_mouse import Mouse_data
from b_load_processed_mat_save_pickle_per_mouse import pickle_dict
from b_load_processed_mat_save_pickle_per_mouse import load_pickleddata



region_data = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/region_based_data_licking_filtered.pickle')
good_regions = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/good_regions.pickle')
phase_dict = {'cond_days':[0,1,2,3,4],
              'deg_days':[5,6,7,8,9],
              'rec_days':[10,11,12,13],
              'ext_days':[14,15,16],
              'finalrec_days':[17,18] }    
#%%
def average_signal_by_trial(multiregion_dict ,region,types,phase,length,num_mouse,session_num,
                            session_of_phase = None,mouse_index = None,rm_bad = True,rm_window = [105,130]):
    key = region
    ave_response_trace = np.ones([session_num,length,num_mouse]) #d1 session, d2 time series, d3 mouse, fill with full trace
    for mouse_id in range(num_mouse):
        for session_of_phase in range(len(phase_dict[phase])):
            x = multiregion_dict[key][types][mouse_id][:,:,phase_dict[phase][session_of_phase]]
            if rm_bad:
                x = remove_bad_trials(x,window = [rm_window[0],rm_window[1]])
            print(key,types,mouse_id)
            ave_response_trace[session_of_phase,:,mouse_id] = np.nanmean(x,axis = 0)
    return ave_response_trace
def remove_bad_trials(mat,window = [105,130],thres = 10):
    mat_window = mat[:,window[0]:window[1]]
    max_val = mat_window.max(axis = 1)
    ind = max_val>=thres
    return mat[ind,:]


length = 180
trialtypes_full = ['go','no_go','UnpredReward','go_omit','background']
average_trace_dict = {}
for phase in phase_dict.keys():
    average_trace_dict[phase] = {}
    for region in good_regions.keys():
        average_trace_dict[phase][region] = {}
        for types in trialtypes_full:
            num_mouse = len(good_regions[region])
            session_num = len(phase_dict[phase])
            if types in ['go','UnpredReward']:
                average_mat = average_signal_by_trial(region_data,region,types,
                                              phase,length,num_mouse,session_num,
                                              rm_bad = False,rm_window = [105,130]) 
            else: 
                average_mat = average_signal_by_trial(region_data,region,types,
                                              phase,length,num_mouse,session_num,
                                              rm_bad = False) 
            # # normalizer always based on go trial
            # go_mat = average_signal_by_trial(region_data,region,'go_omit',
            #                               'cond_days',length,num_mouse,5,
            #                               ) 
            normalizer_peak = np.mean(np.nanmax(average_mat,axis = 1),axis = 0)
            new_mat = average_mat.copy()
            for i in range(session_num):
                for j in range(num_mouse):
                # ax[i].plot(np.nanmean(a[i,:,:],axis =1))
                    new_mat[i,:,j] = average_mat[i,:,j]/normalizer_peak[j]*np.mean(normalizer_peak)
            average_trace_dict[phase][region][types] = new_mat
#%% between regions
#average_trace_dict[phase][region][types]
# 
ttype = 'UnpredReward'
regions = average_trace_dict[phase].keys()
phases = phase_dict.keys()
fig,ax = plt.subplots(len(regions),len(phases) ,figsize = (12,8))
for i,region in enumerate(regions):
    for j,phase in enumerate(phases):
        if ttype == 'go' and phase == 'ext_days':
            ax[i,j].plot(np.nanmean(average_trace_dict[phase][region]['go_omit'][-1,:,:],axis = 1))
        else:
            ax[i,j].plot(np.nanmean(average_trace_dict[phase][region][ttype][-1,:,:],axis = 1))
    
        ax[i,j].spines['top'].set_visible(False)
        ax[i,j].spines['right'].set_visible(False)
        ax[i,j].spines['bottom'].set_visible(False)
        ax[i,j].set_xticks([])
        ax[i,j].spines['left'].set_visible(False)
        ax[i,j].set_ylim([-2,15])
        if j == 0:
            ax[i,j].spines['left'].set_visible(True)
            ax[i,j].set_ylabel(region,fontsize=15)
            ax[i,j].set_yticks(np.array([0,4,8,12]))
            ax[i,j].set_yticklabels(np.array([0,4,8,12]),fontsize=13)
        else:
            ax[i,j].set_yticks([])

        if i == len(regions)-1:
            ax[i,j].spines['bottom'].set_visible(True)
            

            ax[i,j].set_xticks(np.arange(0,180,40))
            
            ax[i,j].set_xticklabels(np.arange(-2,7,2), fontsize=13)
plt.show()    
    
#%%

ttype = 'go_omit'
phase = 'ext_days'
regions = average_trace_dict[phase].keys()

fig,ax = plt.subplots(len(regions),len(phases) ,figsize = (12,8))
for i,region in enumerate(regions):
    for j in range(len(phase_dict[phase])):
        
        ax[i,j].plot(np.nanmean(average_trace_dict[phase][region][ttype][j,:,:],axis = 1))
    
        ax[i,j].spines['top'].set_visible(False)
        ax[i,j].spines['right'].set_visible(False)
        ax[i,j].spines['bottom'].set_visible(False)
        ax[i,j].set_xticks([])
        ax[i,j].spines['left'].set_visible(False)
        ax[i,j].set_ylim([-2,8])
        if j == 0:
            ax[i,j].spines['left'].set_visible(True)
            ax[i,j].set_ylabel(region,fontsize=15)
            ax[i,j].set_yticks(np.array([0,4,8]))
            ax[i,j].set_yticklabels(np.array([0,4,8]),fontsize=13)
        else:
            ax[i,j].set_yticks([])

        if i == len(regions)-1:
            ax[i,j].spines['bottom'].set_visible(True)
            ax[i,j].set_xlabel(f'session {j+1}',fontsize=15)

            ax[i,j].set_xticks(np.arange(0,180,40))
            
            ax[i,j].set_xticklabels(np.arange(-2,7,2), fontsize=13)
plt.show()     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    