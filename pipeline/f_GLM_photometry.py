# -*- coding: utf-8 -*-
"""
Created on Sat Mar 19 17:24:30 2022

@author: qianl
"""
#%% load function and package

from scipy import ndimage as ndi

import pickle
import sys
from skimage import io
import time
import h5py
import os
import math
import pandas as pd
import re
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
from statsmodels.gam.api import GLMGam
import statsmodels.gam.smooth_basis as sb
import statsmodels.api as sm
from patsy import dmatrices
#%% function
def bin_neural_signal(signal,period = [0,180],binframe = 4):
    binned_dff = []
    
    temp_signal = signal[period[0]:period[1]]
    bins = np.arange(45) ## need to adjusted for bins
    for i in bins:
        chunk_dff = signal[i*binframe:(i+1)*binframe]
        mean_in_bin = np.mean(chunk_dff)
        binned_dff.append(mean_in_bin)

    return binned_dff
def get_r2(y,model):
    sst_val = sum(map(lambda x: np.power(x,2),y.values-np.mean(y.values))) 
    sse_val = sum(map(lambda x: np.power(x,2),model.resid_response)) 
    r2 = 1.0 - sse_val/sst_val
    return r2[0]

#%% load data
mouse_name = 'FgDA_01'
gcamp_data = load_pickleddata('D:/PhD/Photometry/DATA/round_20220307/processed/corrected_gcamp_data.pickle')
pickle_data = load_pickleddata('D:/PhD/Photometry/DATA/round_20220307/processed/{}_stats.pickle'.format(mouse_name))
#%%
rel_contrib_DICT = {}
for session in range(len(pickle_data.all_days)):
    if mouse_name == 'FgDA_01':   
        good_sites = ['AMOT','PMOT','ALOT','PLOT','MNacS','LNacS']
    else:
        good_sites = pickle_data.df_bpod_doric[str(session)+'_'+pickle_data.all_days[session]]
    
    num_neuron = len(good_sites)
    rel_contrib_mat = np.zeros([num_neuron,5]) # create a empty matrix to store the relative contribution of each variable
    
    GLM_DICT = {}
    
    for index,site in enumerate(good_sites):
        beh_df = pickle_data.df_bpod_doric[str(session)+'_' + pickle_data.all_days[session]]['dataframe']
        
        binned_licking = []
        binned_dff = []
        diff_binned_licking = []
        trialnum = []
        water = np.zeros([4,45]) 
        go_odor = np.zeros([4,45])
        nogo_odor = np.zeros([4,45])
        
        GLM_df = pd.DataFrame()
        # fluo_signal = gcamp_data[mouse_name][TT]['signal_gcamp']
        for i, row in beh_df.iterrows():
            count_go = 0
            count_nogo = 0
            count_omit = 0
            count_unpred = 0
            if row['Trialtype'] != 'background':
                
                try:
                    binned_licking += list(np.histogram(row['lickings'][0],bins = 45,range = (0,9))[0]/(1/5))
        
                    diff_licking = np.diff(np.histogram(row['lickings'][0],bins = 45,range = (0,9))[0]/(1/5))
                except IndexError:
                    binned_licking += list(np.histogram(row['lickings'],bins = 45,range = (0,9))[0]/(1/5))
        
                    diff_licking = np.diff(np.histogram(row['lickings'],bins = 45,range = (0,9))[0]/(1/5))               
                diff_licking = np.insert(diff_licking,0,0)
                diff_binned_licking += list(diff_licking)
                
                # whole-trial variable
                trialnum += [i]*45
                # events variable
                x = np.linspace(0,9,45)
                knots = sb.get_knots_bsplines(x,df = 4,)#spacing = 'equal')
        
                basis = sb._eval_bspline_basis(x,degree = 3,knots=knots)[0]
                temp_water = np.zeros([4,45])
                temp_go_odor = np.zeros([4,45])
                temp_nogo_odor = np.zeros([4,45])
                if row['Trialtype'] == 'go':
                    
                    signal = list(np.nan_to_num(gcamp_data[mouse_name]['go']['signal_gcamp'][site][count_go,:,session])) # replace nan value with zero and convert to 45 element
                    binned_dff += bin_neural_signal(signal,period = [0,180],binframe = 4)
                    water_event = np.zeros(45)
                    water_event[27] = 1
                    go_odor_event = np.zeros(45)
                    go_odor_event[9:14] = np.ones(5)
                    
                    for i in range(4):
                        conv_water = np.convolve(basis[:,i],water_event,mode = 'full')
                        temp_water[i,:] = conv_water[0:45]
                        conv_go = np.convolve(basis[:,i],go_odor_event,mode = 'full')
                        temp_go_odor[i,:] = conv_go[0:45]
                    count_go += 1
                    
                elif row['Trialtype'] == 'go_omit':
                    signal = list(np.nan_to_num(gcamp_data[mouse_name]['go_omit']['signal_gcamp'][site][count_omit,:,session]))
                    binned_dff += bin_neural_signal(signal,period = [0,180],binframe = 4)
                    go_odor_event = np.zeros(45)
                    go_odor_event[9:14] = np.ones(5)
    
                    for i in range(4):
                        conv_go = np.convolve(basis[:,i],go_odor_event,mode = 'full')
                        temp_go_odor[i,:] = conv_go[0:45]
                    count_omit += 1
                    
                elif row['Trialtype'] == 'no_go':
                    signal = list(np.nan_to_num(gcamp_data[mouse_name]['no_go']['signal_gcamp'][site][count_nogo,:,session]))
                    # water_event = np.array([0,0,0,0,0,0,0])
                    # nogo_odor_event = np.array([0,1,0,0,0,0,0])
                    binned_dff += bin_neural_signal(signal,period = [0,180],binframe = 4)
                    nogo_odor_event = np.zeros(45)
                    nogo_odor_event[9:14] = np.ones(5)
                    for i in range(4):
    
                        conv_nogo = np.convolve(basis[:,i],nogo_odor_event,mode = 'full')
                        temp_nogo_odor[i,:] = conv_nogo[0:45]
                    count_nogo += 1
                elif row['Trialtype'] == 'UnpredReward':
                    signal = list(np.nan_to_num(gcamp_data[mouse_name]['UnpredReward']['signal_gcamp'][site][count_omit,:,session]))
                    # water_event = np.array([0,0,0,0,0,0,0])
                    # nogo_odor_event = np.array([0,1,0,0,0,0,0])
                    binned_dff += bin_neural_signal(signal,period = [0,180],binframe = 4)
                    water_event = np.zeros(45)
                    water_event[27] = 1
                    for i in range(4):
    
                        conv_water = np.convolve(basis[:,i],water_event,mode = 'full')
                        temp_water[i,:] = conv_water[0:45]
                    count_unpred += 1
                water = np.hstack((water,temp_water))
                go_odor = np.hstack((go_odor,temp_go_odor))
                nogo_odor = np.hstack((nogo_odor,temp_nogo_odor))
        
        water = water[:,45:]
        go_odor = go_odor[:,45:]
        nogo_odor = nogo_odor[:,45:]
        GLM_df = pd.DataFrame({'NeuroSignal_dff':binned_dff,'licking':binned_licking,'trialnum':trialnum,'diff_licking':diff_binned_licking,'square_licking':np.power(binned_licking,2)})
        
        for i in range(4):
            GLM_df['water{}'.format(i)] = water[i,:]
        for i in range(4):
            GLM_df['go_odor{}'.format(i)] = go_odor[i,:]    
        for i in range(4):
            GLM_df['nogo_odor{}'.format(i)] = nogo_odor[i,:]    
        GLM_DICT[site] = GLM_df
        
        
            #  1. loop neuron 2. train set and test set (function) 3. fit all variables and partial model, generate R2 4. calculate the relative contribution and build the matr
        # all variables
        formula = 'NeuroSignal_dff ~ licking + diff_licking + square_licking + trialnum + go_odor0 + go_odor1 + go_odor2+go_odor3 + water0+ water1+ water2+ water3+ nogo_odor0+ nogo_odor1+ nogo_odor2+ nogo_odor3 '
        y, X = dmatrices(formula, data=GLM_df, return_type='dataframe')
        glm = sm.GLM(y,X)
        res_o = glm.fit()
        r2= get_r2(y,res_o)
        
        # no licking
        formula = 'NeuroSignal_dff ~ trialnum + go_odor0 + go_odor1 + go_odor2+go_odor3  + water0+ water1+ water2+ water3+ nogo_odor0+ nogo_odor1+ nogo_odor2+ nogo_odor3 '
        y, X = dmatrices(formula, data=GLM_df, return_type='dataframe')
        glm = sm.GLM(y,X)
        res_o_licking = glm.fit()
        r2_licking= get_r2(y,res_o_licking)
        print(r2_licking)
        
        # no go odor
        formula = 'NeuroSignal_dff ~ licking + diff_licking + square_licking + trialnum + water0+ water1+ water2+ water3+  nogo_odor0+ nogo_odor1++ nogo_odor2+ nogo_odor3'
        y, X = dmatrices(formula, data=GLM_df, return_type='dataframe')
        glm = sm.GLM(y,X)
        res_o_go = glm.fit()
        r2_go_odor= get_r2(y,res_o_go)
        
        # no water
        formula = 'NeuroSignal_dff ~ licking+ diff_licking + square_licking + trialnum + go_odor0 + go_odor1 + go_odor2+go_odor3  + nogo_odor0+ nogo_odor1+ nogo_odor2+ nogo_odor3 '
        y, X = dmatrices(formula, data=GLM_df, return_type='dataframe')
        glm = sm.GLM(y,X)
        res_o_water = glm.fit()
        r2_water= get_r2(y,res_o_water)
        
        # no nogo odor
        formula = 'NeuroSignal_dff ~ licking+ diff_licking + square_licking + trialnum + go_odor0 + go_odor1 + go_odor2+go_odor3 + water0+ water1+ water2+ water3'
        y, X = dmatrices(formula, data=GLM_df, return_type='dataframe')
        y, X = dmatrices(formula, data=GLM_df, return_type='dataframe')
        glm = sm.GLM(y,X)
        res_o_nogo = glm.fit()
        r2_nogo_odor= get_r2(y,res_o_nogo)
        
        # no trialnum
        formula = 'NeuroSignal_dff ~ licking+ diff_licking + square_licking + go_odor0 + go_odor1 + go_odor2+go_odor3 + water0+ water1+ water2+ water3+ nogo_odor0+ nogo_odor1+ nogo_odor2+ nogo_odor3 '
        y, X = dmatrices(formula, data=GLM_df, return_type='dataframe')
        glm = sm.GLM(y,X)
        res_o_trialnum = glm.fit()
        r2_trialnum= get_r2(y,res_o_trialnum)
        
        rel_contrib = (r2-r2_go_odor)/((r2-r2_water)+(r2-r2_licking)+(r2-r2_go_odor)+(r2-r2_nogo_odor)+(r2-r2_trialnum))
        print('relative contribution:',rel_contrib)
        rel_contrib_mat[index,0] = rel_contrib
        
        rel_contrib = (r2-r2_nogo_odor)/((r2-r2_water)+(r2-r2_licking)+(r2-r2_go_odor)+(r2-r2_nogo_odor)+(r2-r2_trialnum))
        print('relative contribution:',rel_contrib)
        rel_contrib_mat[index,1] = rel_contrib
        
        rel_contrib = (r2-r2_water)/((r2-r2_water)+(r2-r2_licking)+(r2-r2_go_odor)+(r2-r2_nogo_odor)+(r2-r2_trialnum))
        print('relative contribution:',rel_contrib)
        rel_contrib_mat[index,2] = rel_contrib
        
        rel_contrib = (r2-r2_licking)/((r2-r2_water)+(r2-r2_licking)+(r2-r2_go_odor)+(r2-r2_nogo_odor)+(r2-r2_trialnum))
        print('relative contribution:',rel_contrib)
        rel_contrib_mat[index,3] = rel_contrib
        
        rel_contrib = (r2-r2_trialnum)/((r2-r2_water)+(r2-r2_licking)+(r2-r2_go_odor)+(r2-r2_nogo_odor)+(r2-r2_trialnum))
        print('relative contribution:',rel_contrib)
        rel_contrib_mat[index,4] = rel_contrib
        
    rel_contrib_DICT[str(session)] = rel_contrib_mat


#%%
import seaborn as sns


for session in range(len(pickle_data.all_days)):
    rel_contrib_mat = rel_contrib_DICT[str(session)]
    mean_var = np.mean(rel_contrib_mat,axis = 0)
    std_var = np.std(rel_contrib_mat,axis = 0)
    fig, ax = plt.subplots(figsize = (5,4))
    x = [1,2,3,4,5]
    # ax.bar(x, mean_var, yerr=std_var/2, align='center', alpha=0.5, ecolor='black', capsize=10)
    ax.boxplot(rel_contrib_mat)
    ax.set_ylabel('relatie contribution (%)')
    ax.set_xticks(x)
    ax.set_xticklabels(['go odor','no-go odor','water','licking','trialnum'],rotation = 'vertical')
    ax.set_title('Session {}'.format(session),pad = 20)
    ax.set_ylim([0,1])
    ax.yaxis.grid(True)
    
    # Save the figure and show
    plt.tight_layout()
    # savepath = 'D:/PhD/Photometry/DATA/round_20220307/figures'
    # plt.savefig("{0}/{1}_{2}.png".format(savepath,mouse_name,TT), bbox_inches="tight", dpi = 72)
    # savename = 'D:/PhD/Microscope/Selina/imaging_data/new_figures/{}-{}/relative_contribution_variables_session{}'.format(mouse_name,date,session)
    # plt.savefig(savename+'.png', bbox_inches="tight", dpi = 400,transparent = True)
    
    # plt.savefig(savename+'.svg', bbox_inches="tight", dpi = 400,transparent = True)
    plt.show()
    
    
    
    
    #
    fig, ax = plt.subplots(1,5,figsize = (8,3),sharey = True)
    xlabels = ['go odor','no-go odor','water','licking','trialnum']
    bin_edge = np.linspace(0,1,20)
    yticks = np.arange(0,num_neuron,10)
    for i in range(5):
    
        ax[i].hist(rel_contrib_mat[:,i],bins = bin_edge,alpha=0.5, edgecolor='black')
        
        ax[i].set_xticks([0,1])
        ax[i].set_xticklabels([0,100])
        ax[i].set_yticks(yticks)
        ax[i].set_title(xlabels[i])
        
        
        ax[i].yaxis.grid(True)
    
    # Save the figure and show
    fig.text(0.5, 0.001, 'Relative contribution (%)', ha='center')
    fig.text(0.001, 0.5, 'Neurons', va='center', rotation='vertical')
    plt.tight_layout()
    # savename = 'D:/PhD/Microscope/Selina/imaging_data/new_figures/{}-{}/relative_contribution_variables_histogram_session{}'.format(mouse_name,date,session)
    # plt.savefig(savename+'.png', bbox_inches="tight", dpi = 400,transparent = True)
    
    # plt.savefig(savename+'.svg', bbox_inches="tight", dpi = 400,transparent = True)
    plt.show()
    #
    yticks = np.arange(0,num_neuron,10)
    ax = sns.heatmap(rel_contrib_mat)
    ax.set_ylabel('Neurons')
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    ax.set_xticks(np.array(x)-0.5)
    ax.set_xticklabels(['go odor','no-go odor','water','licking','trialnum'],rotation = 'vertical')
    # savename = 'D:/PhD/Microscope/Selina/imaging_data/new_figures/{}-{}/relative_contribution_variables_heatmap_session{}'.format(mouse_name,date,sessions[session_index])
    # plt.savefig(savename+'.png', bbox_inches="tight", dpi = 400,transparent = True)
    
    # plt.savefig(savename+'.svg', bbox_inches="tight", dpi = 400,transparent = True)
    plt.show()

#%%
contrib_go = []
contrib_nogo = []
contrib_licking = []
contrib_water = []
contrib_trianum = []

contrib_go_std = []
contrib_nogo_std = []
contrib_licking_std = []
contrib_water_std = []
contrib_trianum_std = []

for session in range(len(pickle_data.all_days)):
    rel_contrib_mat = rel_contrib_DICT[str(session)]
    mean_var = np.mean(rel_contrib_mat,axis = 0)   
    std_var = np.std(rel_contrib_mat,axis = 0)
    # adding mean of sites value to variable list
    contrib_go.append(mean_var[0])
    contrib_nogo.append(mean_var[1])
    contrib_licking.append(mean_var[3])
    contrib_water.append(mean_var[2])
    contrib_trianum.append(mean_var[4])
    
    contrib_go_std.append(std_var[0])
    contrib_nogo_std.append(std_var[1])
    contrib_licking_std.append(std_var[3])
    contrib_water_std.append(std_var[2])
    contrib_trianum_std.append(std_var[4])
    if session in [4,9,13,16]:
        contrib_go.append(np.nan)
        contrib_nogo.append(np.nan)
        contrib_licking.append(np.nan)
        contrib_water.append(np.nan)
        contrib_trianum.append(np.nan)  
        
        contrib_go_std.append(np.nan)
        contrib_nogo_std.append(np.nan)
        contrib_licking_std.append(np.nan)
        contrib_water_std.append(np.nan)
        contrib_trianum_std.append(np.nan)
    
#%%

mpl.rcParams['lines.linewidth'] = 1


plt.subplots(1,1,figsize = (11,5))
plt.set_cmap('Set2')
x = np.arange(len(contrib_go))
plt.scatter(x,contrib_go)
plt.errorbar(x,contrib_go,yerr = contrib_go_std,label = 'go') 
plt.scatter(x,contrib_nogo)
plt.errorbar(x,contrib_nogo,yerr = contrib_nogo_std,label = 'no go')   
plt.scatter(x,contrib_water)
plt.errorbar(x,contrib_water,yerr = contrib_water_std,label = 'water')   
plt.scatter(x,contrib_licking)
plt.errorbar(x,contrib_licking,yerr = contrib_licking_std,label = 'licking')  
plt.scatter(x,contrib_trianum) 
plt.errorbar(x,contrib_trianum,yerr = contrib_trianum_std,label = 'motivation')    
plt.legend()
plt.xticks(ticks = x, labels = ['cond1','cond2','cond3','cond4','cond5','','deg1','deg2','deg3','deg4','deg5','','rec1','rec2','rec3','pre-ext','','ext1','ext2','ext3','','frec1','frec2','frec3'],
           rotation = 45)
plt.show()  
    
    