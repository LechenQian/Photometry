# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 19:20:43 2022

@author: qianl


"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
import os
import random
import matplotlib as mpl
import seaborn as sns
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
import matplotlib.cm as cm  
#%%

region_data = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/region_based_data.pickle')
trialtypes = ['go','go_reward','no_go','UnpredReward','UnpredReward_reward','go_omit','background']
good_regions = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/good_regions.pickle')
phase_dict = {'cond_days':[0,1,2,3,4],
              'deg_days':[5,6,7,8,9],
              'rec_days':[10,11,12,13],
              'ext_days':[14,15,16],
              'finalrec_days':[17,18] }    
#%%
 
region = 'LNacS'
types = 'UnpredReward_reward'
phase = 'deg_days'
session_of_phase = 0
mouse_index = 0
length = 40
num_mouse = 5
session_num = 5
def average_signal_by_trial(multiregion_dict ,region,types,phase,length,num_mouse,session_num,
                            session_of_phase = None,mouse_index = None,rm_bad = True,rm_window = [105,130]):
    key = region
    ave_response_trace = np.ones([session_num,length,num_mouse]) #d1 session, d2 time series, d3 mouse, fill with full trace
    for mouse_id in range(num_mouse):
        for session_of_phase in range(len(phase_dict[phase])):
            x = multiregion_dict[key][types][mouse_id][:,:,phase_dict[phase][session_of_phase]]
            if rm_bad:
                x = remove_bad_trials(x,window = [rm_window[0],rm_window[1]])
            ave_response_trace[session_of_phase,:,mouse_id] = np.nanmean(x,axis = 0)
            
    return ave_response_trace
def find_max_response_from_ave_signal(mat):  
    return np.nanmax(mat,axis = 1)


average_mat = average_signal_by_trial(region_data,region,types,
                                      phase,length,num_mouse,session_num,rm_bad = False)   

#%% filter bad trials when there's no water signal
#only keep the trace when the animal actually gets water
def remove_bad_trials(mat,window = [105,130],thres = 10):
    mat_window = mat[:,window[0]:window[1]]
    max_val = mat_window.max(axis = 1)
    ind = max_val>=thres
    return mat[ind,:]
    

# ave_response = find_max_response_from_ave_signal(average_mat)
#%% normalize
normalizer_peak = np.sum(np.nanmax(average_mat,axis = 1),axis = 0)
new_mat = average_mat.copy()
for i in range(5):
    for j in range(5):
    # ax[i].plot(np.nanmean(a[i,:,:],axis =1))
        new_mat[i,:,j] = average_mat[i,:,j]/normalizer_peak[j]#*np.mean(normalizer_peak)



fig,ax = plt.subplots(1,5,sharex=True, sharey=True,figsize = (8,3))

for i in range(5):

    ax[i].plot(np.nanmean(new_mat[i,:,:],axis =1),alpha = 0.8,color = 'k')
    plt.xticks([0,19,39],[0,1,2])
    ax[i].set_xlabel(f'cond{i+1}')
    if i == 0:
        ax[i].set_ylabel('normalized DA signal(Z-score)')        
#%%
fig,ax = plt.subplots(5,5,sharex = True, sharey = True)
for i in range(5):
    for j in range(5):
        ax[i,j].plot(new_mat[i,:,j])
#%% plot last day conditioning and last day degradation
fig = plt.figure()


session_of_phase = 0
mouse_index = 0
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
                                              rm_bad = True,rm_window = [105,130]) 
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



#%% plot

# condition, region LNacS, session 5, go trial



def plot_single_session_average_trace(data, phase,region,types,session_id,CS,US,figsize = (5,4),ylim=[-1,7]):
    aa = np.nanmean(data[phase][region][types][session_id,:,:],axis = 1)
    fig,ax = plt.subplots(figsize = (figsize[0],figsize[1]))
    ax.plot(aa,color = 'k')
    plt.xticks(np.arange(0,180,20),np.arange(-2,7,1))
    plt.xlabel('Time from odor onset(s)')
    plt.ylabel('Signals(Zscores)')
    ax.set_ylim([ylim[0],ylim[1]])
    
    # filled area
    std = np.nanstd(data[phase][region][types][session_id,:,:],axis = 1)
    ax.fill_between(np.arange(0,180,1), aa-std, aa+std,alpha = 0.2,color = 'k')
    ymin, ymax = ax.get_ylim()
    if CS:       
        # vertical lines
        ax.vlines(x=40, ymin=ymin, ymax=ymax, colors='tab:orange', ls='--', lw=2)
    if US:     
        # vertical lines
        ax.vlines(x=110, ymin=ymin, ymax=ymax, colors='tab:blue', ls='--', lw=2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.show()

ymax = 20
plot_single_session_average_trace(data=average_trace_dict, phase = 'cond_days',
                                  region = 'LNacS',types = 'go',session_id = 4,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])

plot_single_session_average_trace(data=average_trace_dict, phase = 'deg_days',
                                  region = 'LNacS',types = 'go',session_id = 4,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])

plot_single_session_average_trace(data=average_trace_dict, phase = 'rec_days',
                                  region = 'LNacS',types = 'go',session_id = 3,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])

plot_single_session_average_trace(data=average_trace_dict, phase = 'ext_days',
                                  region = 'LNacS',types = 'go_omit',session_id = 2,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])

plot_single_session_average_trace(data=average_trace_dict, phase = 'finalrec_days',
                                  region = 'LNacS',types = 'go',session_id = 1,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])



#%%
def plot_multi_session_average_trace(data, phase,region,types,CS=None,US = None,figsize = (4,7),ylim=[-1,7],gap = 7):
     
    fig,ax = plt.subplots(figsize = (figsize[0],figsize[1]))
    i = 0
    sessions = data[phase][region][types].shape[0]
    for session_id in range(sessions):
        aa = np.nanmean(data[phase][region][types][session_id,:,:],axis = 1)
        ax.plot(aa+i,color = cm.cool(session_id/float(sessions)))
        std = np.nanstd(data[phase][region][types][session_id,:,:],axis = 1)
        ax.fill_between(np.arange(0,180,1), aa+i-std, aa+i+std,alpha = 0.2,color = cm.cool(session_id/float(sessions)))
        i += gap
    plt.xticks(np.arange(0,180,20),np.arange(-2,7,1))
    plt.xlabel('Time from odor onset(s)')
    plt.ylabel('Signals(Zscores)')
    # ax.set_ylim([ylim[0],ylim[1]])
    
    # filled area
    
    ymin, ymax = ax.get_ylim()
    if CS:       
        # vertical lines
        ax.vlines(x=40, ymin=ymin, ymax=ymax, colors='tab:orange', ls='--', lw=2)
    if US:     
        # vertical lines
        ax.vlines(x=110, ymin=ymin, ymax=ymax, colors='tab:blue', ls='--', lw=2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([0,5])
    plt.show()
    
plot_multi_session_average_trace(data=average_trace_dict, phase = 'cond_days',
                                  region = 'LNacS',types = 'go',
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap =16)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'deg_days',
                                  region = 'LNacS',types = 'go',
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap = 12)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'rec_days',
                                  region = 'LNacS',types = 'go',
                                  CS=None,US = None,figsize = (4,4.8),
                                  ylim=[-1,7],gap = 14)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'ext_days',
                                  region = 'LNacS',types = 'go_omit',
                                  CS=None,US = None,figsize = (4,3.6),
                                  ylim=[-1,7],gap = 12)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'finalrec_days',
                                  region = 'LNacS',types = 'go',
                                  CS=None,US = None,figsize = (4,2.4),
                                  ylim=[-1,7],gap = 12)

#%% plot multiple trialtypes

def plot_multi_session_multitype_average_trace(data, phase,region,types,CS=None,US = None,figsize = (4,7),ylim=[-1,7],gap = 7):
    type_color = {'go':'orange','go_omit':'yellow','no_go':'green','background':'grey','UnpredReward':'blue'}
    fig,ax = plt.subplots(figsize = (figsize[0],figsize[1]))
    i = 0
    sessions = data[phase][region][types[0]].shape[0]
    for session_id in range(sessions):
        for ttype in types:
            aa = np.nanmean(data[phase][region][ttype][session_id,:,:],axis = 1)
            ax.plot(aa+i,color = type_color[ttype], label = ttype if session_id == 0 else '')
            # std = np.nanstd(data[phase][region][ttype][session_id,:,:],axis = 1)
            # ax.fill_between(np.arange(0,180,1), aa+i-std, aa+i+std,alpha = 0.2,color = type_color[ttype])
        i += gap
    plt.xticks(np.arange(0,180,20),np.arange(-2,7,1))
    plt.xlabel('Time from odor onset(s)')
    plt.ylabel('Signals(Zscores)')
    # ax.set_ylim([ylim[0],ylim[1]])
    
    # filled area
    
    ymin, ymax = ax.get_ylim()
    if CS:       
        # vertical lines
        ax.vlines(x=40, ymin=ymin, ymax=ymax, colors='tab:orange', ls='--', lw=2)
    if US:     
        # vertical lines
        ax.vlines(x=110, ymin=ymin, ymax=ymax, colors='tab:blue', ls='--', lw=2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([0,5])
    plt.legend(bbox_to_anchor=(1.1, 0.5),frameon = False,loc = 'center left')
    plt.show()


plot_multi_session_multitype_average_trace(data=average_trace_dict, phase = 'cond_days',
                                  region = 'LNacS',types = ['go','go_omit','no_go'],
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap = 10)

plot_multi_session_multitype_average_trace(data=average_trace_dict, phase = 'deg_days',
                                  region = 'LNacS',types = ['go','go_omit','no_go','UnpredReward'],
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap = 10)

plot_multi_session_multitype_average_trace(data=average_trace_dict, phase = 'rec_days',
                                  region = 'LNacS',types = ['go','go_omit','no_go'],
                                  CS=None,US = None,figsize = (4,4.8),
                                  ylim=[-1,7],gap = 10)

plot_multi_session_multitype_average_trace(data=average_trace_dict, phase = 'ext_days',
                                  region = 'LNacS',types = ['go_omit','no_go'],
                                  CS=None,US = None,figsize = (4,3.6),
                                  ylim=[-1,7],gap = 10)

plot_multi_session_multitype_average_trace(data=average_trace_dict, phase = 'finalrec_days',
                                  region = 'LNacS',types = ['go','go_omit','no_go'],
                                  CS=None,US = None,figsize = (4,2.4),
                                  ylim=[-1,7],gap = 10)



#%% plotting
# go omit
plot_multi_session_average_trace(data=average_trace_dict, phase = 'cond_days',
                                  region = 'LNacS',types = 'go_omit',
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap =16)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'deg_days',
                                  region = 'LNacS',types = 'go_omit',
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap = 12)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'rec_days',
                                  region = 'LNacS',types = 'go_omit',
                                  CS=None,US = None,figsize = (4,4.8),
                                  ylim=[-1,7],gap = 14)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'ext_days',
                                  region = 'LNacS',types = 'go_omit',
                                  CS=None,US = None,figsize = (4,3.6),
                                  ylim=[-1,7],gap = 12)

plot_multi_session_average_trace(data=average_trace_dict, phase = 'finalrec_days',
                                  region = 'LNacS',types = 'go_omit',
                                  CS=None,US = None,figsize = (4,2.4),
                                  ylim=[-1,7],gap = 12)

#%% plotting
# unpredReward


plot_multi_session_average_trace(data=average_trace_dict, phase = 'deg_days',
                                  region = 'LNacS',types = 'UnpredReward',
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap = 12)




plot_single_session_average_trace(data=average_trace_dict, phase = 'deg_days',
                                  region = 'LNacS',types = 'UnpredReward',session_id = 4,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])



#%% notstacked
def plot_multi_session_average_trace_overlay(data, phase,region,types,CS=None,US = None,figsize = (4,7),ylim=[-1,7],gap = 7):
     
    fig,ax = plt.subplots(figsize = (figsize[0],figsize[1]))
    i = 0
    sessions = data[phase][region][types].shape[0]
    for session_id in range(sessions):
        aa = np.nanmean(data[phase][region][types][session_id,:,:],axis = 1)
        ax.plot(aa+i,color = cm.cool(session_id/float(sessions)))
        # std = np.nanstd(data[phase][region][types][session_id,:,:],axis = 1)
        # ax.fill_between(np.arange(0,180,1), aa+i-std, aa+i+std,alpha = 0.2,color = cm.cool(session_id/float(sessions)))
    plt.xticks(np.arange(0,180,20),np.arange(-2,7,1))
    plt.xlabel('Time from odor onset(s)')
    plt.ylabel('Signals(Zscores)')
    # ax.set_ylim([ylim[0],ylim[1]])
    
    # filled area
    
    ymin, ymax = ax.get_ylim()
    if CS:       
        # vertical lines
        ax.vlines(x=40, ymin=ymin, ymax=ymax, colors='tab:orange', ls='--', lw=2)
    if US:     
        # vertical lines
        ax.vlines(x=110, ymin=ymin, ymax=ymax, colors='tab:blue', ls='--', lw=2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_yticks([0,5])
    plt.show()
        
plot_multi_session_average_trace_overlay(data=average_trace_dict, phase = 'deg_days',
                                  region = 'MNacS',types = 'UnpredReward',
                                  CS=None,US = None,figsize = (4,6),
                                  ylim=[-1,7],gap = 12)