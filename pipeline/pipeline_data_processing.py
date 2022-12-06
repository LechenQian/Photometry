# -*- coding: utf-8 -*-
"""
Created on Sun Oct  2 17:35:34 2022

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
from b_load_processed_mat_save_pickle_per_mouse import Mouse_data
from b_load_processed_mat_save_pickle_per_mouse import pickle_dict
from b_load_processed_mat_save_pickle_per_mouse import load_pickleddata
from sklearn.linear_model import LinearRegression
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

#%%

# mice = ['FgDA_01','FrgD1_01','FrgD2_01','FrgD2_02']
# mice = ['FgDA_02','FgDA_03','FgDA_04','FgDA_05']
mice = ['FgDA_01']
         
         
path = 'D:/PhD/Photometry/DATA/round_20220307'

trialtypes = ['go', 'no_go', 'go_omit', 'background','UnpredReward']

data = {}


# get trials


for mouse_id in mice:
    load_path = os.path.join(path,'processed/{0}_stats.pickle'.format(mouse_id))
    mouse = load_pickleddata(load_path)
    
    #event plot with trials and iscorerct data
    
    # assign two df 
    mouse_trials = mouse.df_bpod_doric
    
    # choose a date
    all_days = mouse.all_days.copy()
    data[mouse_id] = {}
    for trialtype in trialtypes:
        data[mouse_id][trialtype] = {}
        data[mouse_id][trialtype]['signal_gcamp'] = {}
        data[mouse_id][trialtype]['lickings'] =[]
        data[mouse_id][trialtype]['num_trials'] =[]
        
        
        for i,site in enumerate([key for key in mouse_trials[str(0)+'_'+all_days[0]]['corr_ROI'].keys()]):
            data[mouse_id][trialtype]['signal_gcamp'][site] = np.full([160,180,len(all_days)], np.nan)
            
            
            
            for index in range(len(all_days)):
                day = all_days[index] 
                dataframe = mouse_trials[str(index)+'_'+day]['dataframe'].copy()
                is_x_TT = dataframe['Trialtype'] == trialtype #or 'go_omit' # index of go trials
                signal_xTT = dataframe[is_x_TT]
                num_trial = len(signal_xTT[site].values)
                if i == 0 :
                    data[mouse_id][trialtype]['lickings'].append( dataframe[is_x_TT].lickings)
                    data[mouse_id][trialtype]['num_trials'].append(num_trial)
                # movement correction
                # print(trialtype,site,all_days[index])
                try: 
                    concat_func_signal_xTT = np.concatenate(signal_xTT[site].values)
                    concat_isos_signal_xTT = np.concatenate(signal_xTT[site+'_isos'].values)
                    # remove nan value
                    
                    ind_nan = np.isnan(concat_func_signal_xTT) + np.isnan(concat_isos_signal_xTT)
                    concat_func_signal_xTT = concat_func_signal_xTT[~ind_nan]
                    concat_isos_signal_xTT = concat_isos_signal_xTT[~ind_nan]
                    
                    reg = LinearRegression().fit(concat_isos_signal_xTT.reshape(-1,1),concat_func_signal_xTT.reshape(-1,1))
                    beta1 = reg.coef_[0]
                    beta0 = reg.intercept_[0]
                    
                    
                    
                    for trial in range(num_trial):
                        # movement correction
                        signal_subtract = beta0 + beta1*signal_xTT[site+'_isos'].values[trial]
                        signal_corrected = signal_xTT[site].values[trial] - signal_subtract
                        #zscore
                        if np.nanstd(signal_corrected[10:30]) == 0:
                            data[mouse_id][trialtype]['signal_gcamp'][site][trial,:,index] = np.zeros([1,180])
                            
                        elif sum(np.isnan(signal_corrected)) == len(signal_corrected):
                            
                            pass
                        else: 
                            values1 = (signal_corrected - np.nanmean(signal_corrected[10:30]))/np.nanstd(signal_corrected[10:30])
                            data[mouse_id][trialtype]['signal_gcamp'][site][trial,:,index] = values1[0:180]
                except:
                    pass
# adding got_water attributes to go and UnpredReward ttypes

for mouse in mice:
    
    for ttype in ['go','UnpredReward']:
        got_water = []
        
        for ind in range(len(data[mouse][ttype]['num_trials'])):
            if data[mouse][ttype]['num_trials'][ind] == 0:
                got_water.append(None)
                continue
    



        
            licking_data = data[mouse][ttype]['lickings'][ind]
            got_water_session = []
            for val in licking_data.items():
                if val[1].shape[0] == 0:
                    got_water_session.append(0)
                    continue
                
                valid_lick = [x for x in val[1][0] if 5.51<=x <=6]
                if len(valid_lick) == 0:
                    got_water_session.append(0)
                else:
                    got_water_session.append(1)
                    
            print(mouse,ttype,ind,np.mean(got_water_session))
            got_water.append(got_water_session)
            
        data[mouse][ttype]['got_water'] = got_water    
            
# filter signal by licking
for mouse in mice:
    
    
    for ttype in ['go','UnpredReward']:
        
        data[mouse][ttype]['signal_gcamp_lick_filtered'] = {}
        for site in data[mouse][ttype]['signal_gcamp'].keys():
            
            session_num = data[mouse][ttype]['signal_gcamp'][site].shape[2]
            data[mouse][ttype]['signal_gcamp_lick_filtered'][site] = np.full([60,180,session_num], np.nan)    
            
            
            for ind in range(session_num):
                
                if data[mouse][ttype]['got_water'][ind] is not None:
                    i = 0
                    for j,val in enumerate(data[mouse][ttype]['got_water'][ind]):
                        if val == 1:
                            data[mouse][ttype]['signal_gcamp_lick_filtered'][site][i,:,ind] = data[mouse][ttype]['signal_gcamp'][site][j,:,ind]
                            i +=1
                
    


save_path = path+'/processed'
filename = 'corrected_gcamp_data_add_licking'
pickle_dict(data,save_path,filename)  



#%%% step 2: update combined data with new batch
data1 = load_pickleddata('D:/PhD/Photometry/DATA/round_20220307/processed/corrected_gcamp_data_add_licking.pickle')
# the following is the one to be added
data2 = load_pickleddata('D:/PhD/Photometry/DATA/round_20220504/processed/corrected_gcamp_data_add_licking.pickle')

data1.update(data2)
#%% step 2.1: save updated dict
save_path = 'D:/PhD/Photometry/results/corrected_DA_gcamp_signals'
filename = 'combined_corrected_DA_gcamp_data_add_licking'
pickle_dict(data1,save_path,filename) 

#%% step 3 change the data structure to brain region based, align the unpredicted reward response and save
# load data
data = load_pickleddata('D:/PhD/Photometry/results/corrected_DA_gcamp_signals/combined_corrected_DA_gcamp_data_add_licking.pickle')
# dict for good regions
good_regions = {'AMOT':['FgDA_01','FgDA_03','FgDA_04'],
                'PMOT':['FgDA_01','FgDA_02','FgDA_04','FgDA_05'],
                'ALOT':['FgDA_03','FgDA_04','FgDA_05'],
                'PLOT':['FgDA_01','FgDA_02','FgDA_04','FgDA_05'],
                'MNacS':['FgDA_01','FgDA_02','FgDA_03','FgDA_04'],
                'LNacS':['FgDA_01','FgDA_02','FgDA_03','FgDA_04','FgDA_05']}
# load_pickleddata('D:/PhD/Photometry/DATA/round_20220504/processed/good_regions.pickle')
#%%
# trialtypes
ttypes= ['go','no_go','UnpredReward','go_omit','background']
multiregion_dict = {}
for key, vals in good_regions.items():
    multiregion_dict[key] = multiregion_dict.get(key,{})
    for val in vals:
        for types in ttypes:
            multiregion_dict[key][types] = multiregion_dict[key].get(types,[])
            if types in ['go','UnpredReward']:
                multiregion_dict[key][types].append(data[val][types]['signal_gcamp_lick_filtered'][key][:,:,:19])    
            else:
                multiregion_dict[key][types].append(data[val][types]['signal_gcamp'][key][:,:,:19])    
    
    

    
save_path = 'D:/PhD/Photometry/results/corrected_DA_gcamp_signals'
filename = 'good_regions'
pickle_dict(good_regions,save_path,filename)  
filename = 'region_based_data_licking_filtered'
pickle_dict(multiregion_dict,save_path,filename)     
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    