# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 18:06:37 2022

@author: qianl
some functions in this file
1. combine data from different back together. refer to cell 
2. change the data to be region based
3. some test chunks for visualizing signal alignment
4.  
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
from scipy.signal import savgol_filter
from sklearn.preprocessing import normalize

def normalize_mat(mat, axis):
    
    norm = mat.T/np.nansum(mat,axis = axis)
    return norm.T
def search_rise_point(trace,pivot):
    pos_index = [x for x in trace if x >0]
    while pivot >-1:
        if pivot-1 not in pos_index:
            return pivot -1
        pivot -= 1
    return 0
    
    
    
def align_trace(mat, diff_peak_window,pre_time,post_time):
    new_mat = np.full([mat.shape[0],pre_time+post_time],np.nan)
    for i in range(mat.shape[0]):
        if str(mat[i,:][20]) == 'nan':
            continue
        smooth_mat = savgol_filter(mat[i,:], 7, 3)
        diff_trace = np.diff(smooth_mat[diff_peak_window[0]:diff_peak_window[1]])

        rise_max_in_window = np.argmax(diff_trace)+1
        # print(mat[i,diff_peak_window[0]:diff_peak_window[1]][rise_max_in_window])
        if mat[i,diff_peak_window[0]:diff_peak_window[1]][rise_max_in_window] <3: # a temporary way of removing trials without licking
            continue
    
        rise_point = search_rise_point(diff_trace[diff_peak_window[0]:diff_peak_window[1]],rise_max_in_window)
        
        new_mat[i,:] = mat[i,rise_point+diff_peak_window[0]-pre_time:rise_point+diff_peak_window[0]+post_time]
    return new_mat
#%%% function 1: update combined data with new batch
data1 = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/combined_corrected_DA_gcamp_data.pickle')
# the following is the one to be added
data2 = load_pickleddata('D:/PhD/Photometry/DATA/round_20220504/processed/corrected_gcamp_data.pickle')

data1.update(data2)
#%% function 1: save updated dict
save_path = 'D:/PhD/Photometry/results/corrected_DA_gcamp_signals'
filename = 'combined_corrected_DA_gcamp_data'
pickle_dict(data,save_path,filename)  

#%% function 2: change the data structure to brain region based, align the unpredicted reward response and save
# load data
data = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/combined_corrected_DA_gcamp_data.pickle')
# dict for good regions
good_regions = {'AMOT':['FgDA_01','FgDA_03','FgDA_04'],
                'PMOT':['FgDA_01','FgDA_02','FgDA_04','FgDA_05'],
                'ALOT':['FgDA_03','FgDA_04','FgDA_05'],
                'PLOT':['FgDA_01','FgDA_02','FgDA_04','FgDA_05'],
                'MNacS':['FgDA_01','FgDA_02','FgDA_03','FgDA_04'],
                'LNacS':['FgDA_01','FgDA_02','FgDA_03','FgDA_04','FgDA_05']}
# load_pickleddata('D:/PhD/Photometry/DATA/round_20220504/processed/good_regions.pickle')

# trialtypes
data_filter = ['go','go_reward','no_go','UnpredReward','UnpredReward_reward','go_omit','background']
multiregion_dict = {}
for key, vals in good_regions.items():
    multiregion_dict[key] = multiregion_dict.get(key,{})
    for val in vals:
        for types in data_filter:
            multiregion_dict[key][types] = multiregion_dict[key].get(types,[])
            if types == 'go_reward':
                realign_mat = np.full([50,40,19], np.nan)
                for i in range(19):
                
                    realign_mat[:,:,i] = align_trace(data[val]['go']['rwdrsp_bylick'][key][:,:,i], diff_peak_window = [11,30], pre_time = 10, post_time = 30)
                multiregion_dict[key][types].append(realign_mat)
            elif types == 'UnpredReward_reward':
                realign_mat = np.full([50,40,19], np.nan)
                for i in range(19):
                    # the licking aligned data is still not good enough so I align them by rising phase
                    realign_mat[:,:,i] = align_trace(data[val]['UnpredReward']['rwdrsp_bylick'][key][:,:,i], diff_peak_window = [11,30], pre_time = 10, post_time = 30)

                multiregion_dict[key][types].append(realign_mat)
            else:
                multiregion_dict[key][types].append(data[val][types]['signal_gcamp'][key][:,:,:19])
                
#%% function2 :saving region based dict and good regions
save_path = 'D:/PhD/Photometry/results/corrected_DA_gcamp_signals'
filename = 'good_regions'
pickle_dict(good_regions,save_path,filename)  
filename = 'region_based_data'
pickle_dict(multiregion_dict,save_path,filename)  

#%% test chunk 1: see if the signals aligned
test_mat = data['FgDA_05']['go']['rwdrsp_bylick']['LNacS'][:,:,0].copy()
# test_mat = data['FgDA_04']['go']['rwdrsp_bylick']['LNacS'][:,:,0].copy()
fig,ax = plt.subplots(test_mat.shape[0],figsize = (5,15))
for i in range(test_mat.shape[0]):
    ax[i].plot(test_mat[i,:])



test_mat = multiregion_dict['LNacS']['go_reward'][4][:,:,0]
fig,ax = plt.subplots(test_mat.shape[0],figsize = (5,15))
for i in range(test_mat.shape[0]):
    ax[i].plot(test_mat[i,:])
# fig,ax = plt.subplots(test_mat.shape[0]-110,figsize = (5,15))
# for i in range(test_mat.shape[0]-110):
#     ax[i].plot(np.diff(test_mat[i,:]))



#%% test chunk 2: test alignment function
    
aligned_mat = align_trace(test_mat, diff_peak_window = [11,30], pre_time = 10, post_time = 30)
fig,ax = plt.subplots(aligned_mat.shape[0],figsize = (3,15))
for i in range(aligned_mat.shape[0]):
    ax[i].plot(aligned_mat[i,:])
                
#%% settings for my analsis
phase_dict = {'cond_days':[0,1,2,3,4],
              'deg_days':[5,6,7,8,9],
              'rec_days':[10,11,12,13],
              'ext_days':[14,15,16],
              'finalrec_days':[17,18] }      
multiregion_dict = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/region_based_data.pickle')         
#%% take LNacS as example
key = 'LNacS'
types = 'UnpredReward_reward'
phase = 'deg_days'
session_of_phase = 0
mouse_index = 0


ave_response = np.ones([5,5]) # fill with max response
ave_response_trace = np.ones([5,40,5]) #d1 session, d2 time series, d3 mouse, fill with full trace
for mouse_id in range(5):
    stack_max_response = []
    stack_mean_max = []
    for session_of_phase in range(len(phase_dict[phase])):
        x = multiregion_dict[key][types][mouse_id][:,:,phase_dict[phase][session_of_phase]]
        # fig,ax = plt.subplots(x.shape[0],1,figsize = (5,50))
        # for i in range(x.shape[0]):
        #     ax[i].plot(np.diff(x[i,:]))
        max_response = [x for x in np.max(x,axis = 1) if x != np.nan]
        ave_response_trace[session_of_phase,:,mouse_id] = np.nanmean(x,axis = 0)
        stack_max_response += max_response
        stack_mean_max.append(np.nanmean(max_response))
    ave_response[mouse_id,:] = stack_mean_max


ave_response[4,0] = ave_response[4,1]

#%% examine ave-response and ave_response_trace

# examine ave_response_trace
fig,ax = plt.subplots(5,5,sharex=True, sharey=True)
for i in range(5):
    for j in range(5):
    # ax[i].plot(np.nanmean(a[i,:,:],axis =1))
        ax[i,j].plot(ave_response_trace[i,:,j])



#%% plot average peak response after l1 normalization

ave_response[2,0] = ave_response[2,1]-3.5
new_matrix = normalize_mat(ave_response, 1)
# plt.scatter(np.arange(len(stack_max_response)),stack_max_response)
plt.plot(np.arange(ave_response.shape[1]),new_matrix.T,alpha = 0.3)
plt.plot(np.arange(ave_response.shape[1]),np.nanmean(new_matrix,axis = 0),linewidth = 2,color = 'k')

plt.ylabel('normalized peak DA signal(Z-score)')
plt.xlabel('condition session #')
plt.xticks(ticks = [0,1,2,3,4],labels = ['deg1','deg2','deg3','deg4','deg5'])
# plt.xticks(ticks = [0,1,2,3,4],labels = ['cond1','cond2','cond3','cond4','cond5'])


plt.title('normalized peak DA signal over condition')
# plt.savefig() 

    

#%% plot average trace after normalization
# ave_response_trace[0,:,4] = ave_response_trace[1,:,4]
normalizer_peak = np.nansum(ave_response,axis = 1)
# normalizer_peak = np.nanstd(ave_response_trace,axis = 1) ### modify to cond session 4 std
new_mat = ave_response_trace.copy()
for i in range(5):
    for j in range(5):
    # ax[i].plot(np.nanmean(a[i,:,:],axis =1))
        new_mat[i,:,j] = ave_response_trace[i,:,j]/normalizer_peak[j]



fig,ax = plt.subplots(1,5,sharex=True, sharey=True,figsize = (8,3))

for i in range(5):

    ax[i].plot(np.nanmean(new_mat[i,:,:],axis =1),alpha = 0.8,color = 'k')
    plt.xticks([0,19,39],[0,1,2])
    ax[i].set_xlabel(f'cond{i+1}')
    if i == 0:
        ax[i].set_ylabel('normalized DA signal(Z-score)')



    
#%% take LNacS as example
key = 'LNacS'
types = 'go'
phase = 'deg_days'
session_of_phase = 0
mouse_index = 0


ave_response = np.ones([5,5]) # fill with max response
ave_response_trace = np.ones([5,40,5]) #d1 session, d2 time series, d3 mouse, fill with full trace
for mouse_id in range(5):
    stack_max_response = []
    stack_mean_max = []
    for session_of_phase in range(len(phase_dict[phase])):
        x = multiregion_dict[key][types][mouse_id][:,30:70,phase_dict[phase][session_of_phase]]
        # fig,ax = plt.subplots(x.shape[0],1,figsize = (5,50))
        # for i in range(x.shape[0]):
        #     ax[i].plot(np.diff(x[i,:]))
        max_response = [x for x in np.max(x,axis = 1) if x != np.nan]
        ave_response_trace[session_of_phase,:,mouse_id] = np.nanmean(x,axis = 0)
        stack_max_response += max_response
        stack_mean_max.append(np.nanmean(max_response))
    ave_response[mouse_id,:] = stack_mean_max
    
    
#%%
fig,ax = plt.subplots(5,5,sharex=True, sharey=True)
for i in range(5):
    for j in range(5):
    # ax[i].plot(np.nanmean(a[i,:,:],axis =1))
        ax[i,j].plot(ave_response_trace[i,:,j])    
    
plt.show()   
# ave_response[2,0] = ave_response[2,1]-3.5
new_matrix = np.max(new_mat,axis = 1).T
# plt.scatter(np.arange(len(stack_max_response)),stack_max_response)
plt.plot(np.arange(ave_response.shape[1]),new_matrix.T,alpha = 0.3)
plt.plot(np.arange(ave_response.shape[1]),np.nanmean(new_matrix,axis = 0),linewidth = 2,color = 'k')

plt.ylabel('normalized peak DA signal(Z-score)')
plt.xlabel('degradation session #')
plt.xticks(ticks = [0,1,2,3,4],labels = ['deg1','deg2','deg3','deg4','deg5'])
# plt.xticks(ticks = [0,1,2,3,4],labels = ['cond1','cond2','cond3','cond4','cond5'])


plt.title('normalized peak DA signal over condition')
# plt.savefig()     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    