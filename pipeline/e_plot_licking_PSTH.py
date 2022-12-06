# -*- coding: utf-8 -*-
"""
Created on Sun Mar 20 23:10:09 2022

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
import scipy.stats as stats
from scipy import signal
from statsmodels.gam.api import GLMGam
import statsmodels.gam.smooth_basis as sb
import statsmodels.api as sm

#%%

# mice = ['FgDA_01','FrgD1_01','FrgD2_01','FrgD2_02']
mice = ['FgDA_02','FgDA_03','FgDA_04','FgDA_05','FgDA_01']
path = 'D:/PhD/Photometry/results/pickles/DA'

trialtypes = ['go', 'no_go', 'go_omit', 'background','UnpredReward']

licking_data = {}


# get trials

sample_points = 45
for mouse_id in mice:
    load_path = os.path.join(path,'{0}_stats.pickle'.format(mouse_id))
    mouse = load_pickleddata(load_path)
    
    #event plot with trials and iscorerct data
    
    # assign two df 
    mouse_trials = mouse.df_bpod_doric
    
    # choose a date
    all_days = mouse.all_days.copy()
    licking_data[mouse_id] = {}
    for trialtype in trialtypes:
        
        
        
    
    
        licking_data[mouse_id][trialtype] = np.full([160,sample_points,len(all_days)], np.nan)
        
        
        
        for index in range(len(all_days)):
            day = all_days[index] 
            dataframe = mouse_trials[str(index)+'_'+day]['dataframe'].copy()
            is_x_TT = dataframe['Trialtype'] == trialtype #or 'go_omit' # index of go trials
            data_xTT = dataframe[is_x_TT]
            num_trial = np.sum(is_x_TT)

            for trial in range(num_trial):
                lickings = data_xTT['lickings'].values[trial]
                if len(lickings) == 0:
                    binned_licks = np.histogram(lickings,bins = sample_points,range = (0,9))[0]/(9/sample_points)
                    
                else:
                    lickings = [i for i in lickings[0] if i <9]
                    binned_licks = np.histogram(lickings,bins = sample_points,range = (0,9))[0]/(9/sample_points)

                
                licking_data[mouse_id][trialtype][trial,:,index] = binned_licks


   #%%         
save_path = 'D:/PhD/Photometry/results/lickings'
filename = 'licking_data'
pickle_dict(licking_data,save_path,filename)           


#%%
import matplotlib.cm as cm
import matplotlib as mpl
rc = {"axes.spines.left" : False,
      "axes.spines.right" : False,
      "axes.spines.bottom" : False,
      "axes.spines.top" : False,
}
plt.rcParams.update(rc)

mpl.rcParams['lines.linewidth'] = 4

cond_days = [0,1,2,3,4]
deg_days = [5,6,7,8,9]
rec_days = [10,11,12,13]
ext_days = [14,15,16]
finalrec_days = [17,18]

sample_points = 45
water_point = int(sample_points/9*5.5)
odor_on_point = int(sample_points/9*2)
odor_off_point = int(sample_points/9*3)

gap = 15

x = np.linspace(0,1,10)
knots = sb.get_knots_bsplines(x,df = 7)#spacing = 'equal')

basis = sb._eval_bspline_basis(x,degree = 3,knots=knots)[0][:,1]

# mice = ['FgDA_01','FrgD1_01','FrgD2_01','FrgD2_02']
mice = ['FgDA_02','FgDA_03','FgDA_04','FgDA_05','FgDA_01']
trialtypes = ['go', 'no_go', 'go_omit', 'UnpredReward','background']

for mouse_name in mice:
    for TT in trialtypes:
        if TT == 'go':
            water_point = int(sample_points/9*5.5)
            odor_on_point = int(sample_points/9*2)
            odor_off_point = int(sample_points/9*3)
        elif TT in ['go_omit','no_go']:
            odor_on_point = int(sample_points/9*2)
            odor_off_point = int(sample_points/9*3)
            water_point = np.nan
        elif TT == 'UnpredReward':
            water_point = int(sample_points/9*5.5)
            odor_on_point = np.nan
            odor_off_point = np.nan
        else:
            water_point = np.nan
            odor_on_point = np.nan
            odor_off_point = np.nan
                    
        
        fig,ax = plt.subplots(1,5,sharex = True, sharey = True,figsize = (12,7))
        
        plt.setp(ax, xticks=np.arange(0,sample_points,int(sample_points/9)), xticklabels=np.arange(0,9,1),
                yticks=np.arange(0,15,10))
        
        
        across_day_licking = np.nanmean(licking_data[mouse_name][TT],axis = 0)
        
        if TT =='go':
            across_day_licking_omit = np.nanmean(licking_data[mouse_name]['go_omit'],axis = 0)
        else:
            across_day_licking_omit = np.nanmean(licking_data[mouse_name][TT],axis = 0)
            
        for index,i in enumerate(cond_days):
            licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
            ax[0].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
            ax[0].legend(frameon=False)
            ax[0].set_ylim(-2,gap*5+20)
            ax[0].set_title('conditioning',fontsize = 20)
            ax[0].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[0].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[0].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
            ax[0].set_ylabel('licking rate (/s)',fontsize = 35)
            
        for index,i in enumerate(deg_days):
            licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
            ax[1].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
            ax[1].legend(frameon=False)
            ax[1].set_title('degradation',fontsize = 20)
            ax[1].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[1].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[1].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
        for index,i in enumerate(rec_days):
            licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
            ax[2].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
            ax[2].legend(frameon=False)
            ax[2].set_title('recovery',fontsize = 20)
            ax[2].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[2].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[2].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
        for index,i in enumerate(ext_days):
            licking_conv = np.convolve(basis,across_day_licking_omit[:,i],mode = 'full')[0:sample_points]
            ax[3].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
            ax[3].legend(frameon=False)
            ax[3].set_title('extinction',fontsize = 20)
            ax[3].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[3].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[3].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
        for index,i in enumerate(finalrec_days):
            licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
            ax[4].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
            ax[4].legend(frameon=False)
            ax[4].set_title('final recovery',fontsize = 20)
            ax[4].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[4].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
            ax[4].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
    
        
        plt.suptitle('trial type{}'.format(TT))
        savepath = 'D:/PhD/Photometry/results/plots/lickings'
        plt.savefig("{0}/{1}_{2}_licking_{3}.png".format(savepath,mouse_name,TT,date.today()), bbox_inches="tight", dpi = 72)
        plt.show()

#%% averaged licking plot






# mice = ['FgDA_01','FrgD1_01','FrgD2_01','FrgD2_02']
mice = ['FgDA_02','FgDA_03','FgDA_04','FgDA_05','FgDA_01']
trialtypes = ['go', 'no_go', 'go_omit', 'UnpredReward','background']


mean_licking_mat = np.zeros([len(mice),sample_points,19])

mean_licking_bymouse_dict = {'go':np.zeros([len(mice),sample_points,19]), 'no_go':np.zeros([len(mice),sample_points,19]), 
                            'go_omit':np.zeros([len(mice),sample_points,19]), 
                             'UnpredReward':np.zeros([len(mice),sample_points,19]),'background':np.zeros([len(mice),sample_points,19])}

for i,mouse_name in enumerate(mice):
    for TT in trialtypes:
        lick_trace = licking_data[mouse_name][TT]
        mean_lick_trace = np.nanmean(lick_trace[:,:,:19],axis = 0)
        mean_licking_bymouse_dict[TT][i,:,:] = mean_lick_trace
        


#%% plot averaged licking traces

import matplotlib.cm as cm
import matplotlib as mpl
rc = {"axes.spines.left" : False,
      "axes.spines.right" : False,
      "axes.spines.bottom" : False,
      "axes.spines.top" : False,
}
plt.rcParams.update(rc)

mpl.rcParams['lines.linewidth'] = 4

cond_days = [0,1,2,3,4]
deg_days = [5,6,7,8,9]
rec_days = [10,11,12,13]
ext_days = [14,15,16]
finalrec_days = [17,18]

sample_points = 45
water_point = int(sample_points/9*5.5)
odor_on_point = int(sample_points/9*2)
odor_off_point = int(sample_points/9*3)

gap = 15

x = np.linspace(0,1,10)
knots = sb.get_knots_bsplines(x,df = 7)#spacing = 'equal')

basis = sb._eval_bspline_basis(x,degree = 3,knots=knots)[0][:,1]

trialtypes = ['go', 'no_go', 'go_omit', 'UnpredReward','background']


for TT in trialtypes:
    if TT == 'go':
        water_point = int(sample_points/9*5.5)
        odor_on_point = int(sample_points/9*2)
        odor_off_point = int(sample_points/9*3)
    elif TT in ['go_omit','no_go']:
        odor_on_point = int(sample_points/9*2)
        odor_off_point = int(sample_points/9*3)
        water_point = np.nan
    elif TT == 'UnpredReward':
        water_point = int(sample_points/9*5.5)
        odor_on_point = np.nan
        odor_off_point = np.nan
    else:
        water_point = np.nan
        odor_on_point = np.nan
        odor_off_point = np.nan
                
    
    fig,ax = plt.subplots(1,5,sharex = True, sharey = True,figsize = (12,7))
    
    plt.setp(ax, xticks=np.arange(0,sample_points,int(sample_points/9)), xticklabels=np.arange(0,9,1),
            yticks=np.arange(0,15,10))
    
    
    across_day_licking = np.nanmean(mean_licking_bymouse_dict[TT],axis = 0)
    
    if TT =='go':
        across_day_licking_omit = np.nanmean(mean_licking_bymouse_dict['go_omit'],axis = 0)
    else:
        across_day_licking_omit = np.nanmean(mean_licking_bymouse_dict[TT],axis = 0)
            
    for index,i in enumerate(cond_days):
        licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
        ax[0].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
        ax[0].legend(frameon=False)
        ax[0].set_ylim(-2,gap*5+20)
        ax[0].set_title('conditioning',fontsize = 20)
        ax[0].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[0].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[0].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
        ax[0].set_ylabel('licking rate (/s)',fontsize = 35)
        
    for index,i in enumerate(deg_days):
        licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
        ax[1].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
        ax[1].legend(frameon=False)
        ax[1].set_title('degradation',fontsize = 20)
        ax[1].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[1].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[1].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
    for index,i in enumerate(rec_days):
        licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
        ax[2].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
        ax[2].legend(frameon=False)
        ax[2].set_title('recovery',fontsize = 20)
        ax[2].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[2].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[2].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
    for index,i in enumerate(ext_days):
        licking_conv = np.convolve(basis,across_day_licking_omit[:,i],mode = 'full')[0:sample_points]
        ax[3].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
        ax[3].legend(frameon=False)
        ax[3].set_title('extinction',fontsize = 20)
        ax[3].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[3].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[3].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
    for index,i in enumerate(finalrec_days):
        licking_conv = np.convolve(basis,across_day_licking[:,i],mode = 'full')[0:sample_points]
        ax[4].plot(licking_conv+gap*index,color =cm.spring(index/float(5)),label = all_days[i])
        ax[4].legend(frameon=False)
        ax[4].set_title('final recovery',fontsize = 20)
        ax[4].axvline(x = odor_on_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[4].axvline(x = odor_off_point,linewidth=0.2, color=(0, 0, 0, 0.75))
        ax[4].axvline(x = water_point,linewidth=0.5, color=(0.17578125, 0.67578125, 0.83203125,1))
    plt.suptitle('trialtype {}'.format(TT))
    savepath = 'D:/PhD/Photometry/results/plots/lickings'
    plt.savefig("{0}/{1}_mean_licking_{2}.png".format(savepath,TT,date.today()), bbox_inches="tight", dpi = 72)
    plt.show()
    
    
#%% plot single date average licking trace

fig = plt.figure()


session_of_phase = 0

length = 45
trialtypes_full = ['go','no_go','UnpredReward','go_omit','background']

def plot_single_session_average_licking(data, types,session_id,CS,US,figsize = (5,4),ylim=[-1,7]):
    aa = np.nanmean(data[types][:,:,session_id],axis = 0)
    fig,ax = plt.subplots(figsize = (figsize[0],figsize[1]))
    ax.plot(aa,color = 'k')
    plt.xticks(np.arange(0,45,5),np.arange(-2,7,1))
    plt.xlabel('Time from odor onset(s)')
    plt.ylabel('licking rate (/s)')
    ax.set_ylim([ylim[0],ylim[1]])
    
    # filled area
    std = np.nanstd(data[types][:,:,session_id],axis = 0)
    ax.fill_between(np.arange(0,45,1), aa-std, aa+std,alpha = 0.2,color = 'k')
    ymin, ymax = ax.get_ylim()
    if CS:       
        # vertical lines
        ax.vlines(x=10, ymin=ymin, ymax=ymax, colors='tab:orange', ls='--', lw=2)
    if US:     
        # vertical lines
        ax.vlines(x=27.5, ymin=ymin, ymax=ymax, colors='tab:blue', ls='--', lw=2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    savepath = 'D:/PhD/Photometry/results/plots/lickings'
    plt.savefig("{0}/{1}_{2}_mean_licking_{3}.png".format(savepath,TT,session_id,date.today()), bbox_inches="tight", dpi = 200)
    plt.show()
    
ymax = 20
TT = 'go_omit'
plot_single_session_average_licking(data = mean_licking_bymouse_dict, types = TT ,session_id = 4,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])  
plot_single_session_average_licking(data = mean_licking_bymouse_dict, types = TT ,session_id = 5,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])  
plot_single_session_average_licking(data = mean_licking_bymouse_dict, types = TT ,session_id = 9,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])  
plot_single_session_average_licking(data = mean_licking_bymouse_dict, types = TT ,session_id = 13,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])  
plot_single_session_average_licking(data = mean_licking_bymouse_dict, types = TT ,session_id = 14,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])  
plot_single_session_average_licking(data = mean_licking_bymouse_dict, types = TT ,session_id = 16,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])  
plot_single_session_average_licking(data = mean_licking_bymouse_dict, types = TT ,session_id = 18,
                                  CS = True,US = True,
                                  figsize = (5,4),ylim=[-1,ymax])  



#%% quantify anticipatory licking
antip_window = [16,25]
fig,ax = plt.subplots(figsize = (8,3))
TT = 'go'

plt.plot(np.arange(0,5),np.nanmean(np.nanmean(mean_licking_bymouse_dict[TT][:,antip_window[0]:antip_window[1],0:5],axis = 1),axis = 0),
         'o-',c= '#A9B9DF',markersize = 8)
plt.plot(np.arange(5,10),np.nanmean(np.nanmean(mean_licking_bymouse_dict[TT][:,antip_window[0]:antip_window[1],5:10],axis = 1),axis = 0),
         'o-',c = '#FFC5CE',markersize = 8)
plt.plot(np.arange(10,14),np.nanmean(np.nanmean(mean_licking_bymouse_dict[TT][:,antip_window[0]:antip_window[1],10:14],axis = 1),axis = 0),
         'o-',c='#A9B9DF',markersize = 8)
if TT == 'go':
    plt.plot(np.arange(14,17),np.nanmean(np.nanmean(mean_licking_bymouse_dict['go_omit'][:,antip_window[0]:antip_window[1],14:17],axis = 1),axis = 0),
         'o-',c='#F6D7C3',markersize = 8)
else:
    plt.plot(np.arange(14,17),np.nanmean(np.nanmean(mean_licking_bymouse_dict[TT][:,antip_window[0]:antip_window[1],14:17],axis = 1),axis = 0),
         'o-',c='#F6D7C3',markersize = 8)
    
plt.plot(np.arange(17,19),np.nanmean(np.nanmean(mean_licking_bymouse_dict[TT][:,antip_window[0]:antip_window[1],17:19],axis = 1),axis = 0),
         'o-',c='#A9B9DF',markersize = 8)

plt.xticks(np.arange(19))
plt.ylabel('mean licking rate (/s)')
plt.xlabel('# sessions')
ax.spines['left'].set_visible(True)
ax.spines['bottom'].set_visible(True)
savepath = 'D:/PhD/Photometry/results/plots/lickings'
plt.savefig("{0}/{1}_mean_licking_all_phases_{2}.png".format(savepath,TT,date.today()), 
            bbox_inches="tight", dpi = 72)
plt.show()

#%% load lickingdata
licking_data = load_pickleddata('D:/PhD/Photometry/results/lickings/licking_data.pickle')


lickings1,lickings2,lickings3,lickings4,lickings5 = [],[],[],[],[]
trial_nums1,trial_nums2,trial_nums3,trial_nums4,trial_nums5 = [],[],[],[],[]
for i in range(5,10):
    mat1 = licking_data['FgDA_01']['background'][:,:,i]
    mat2 = licking_data['FgDA_02']['background'][:,:,i]
    mat3 = licking_data['FgDA_03']['background'][:,:,i]
    mat4 = licking_data['FgDA_04']['background'][:,:,i]
    mat5 = licking_data['FgDA_05']['background'][:,:,i]
    ind_val = [ not x for x in np.isnan(mat1[:,0])]
    trial_nums1.append(np.sum(ind_val))
    lickings1+=list(np.mean(mat1[ind_val,10:20],axis = 1))
    
    ind_val = [ not x for x in np.isnan(mat2[:,0])]
    trial_nums2.append(np.sum(ind_val))
    lickings2+=list(np.mean(mat2[ind_val,10:20],axis = 1))
    
    ind_val = [ not x for x in np.isnan(mat3[:,0])]
    trial_nums3.append(np.sum(ind_val))
    lickings3+=list(np.mean(mat3[ind_val,10:20],axis = 1))
    
    ind_val = [ not x for x in np.isnan(mat4[:,0])]
    trial_nums4.append(np.sum(ind_val))
    lickings4+=list(np.mean(mat4[ind_val,10:20],axis = 1))
    
    ind_val = [ not x for x in np.isnan(mat5[:,0])]
    trial_nums5.append(np.sum(ind_val))
    lickings5+=list(np.mean(mat5[ind_val,10:20],axis = 1))



fig,ax = plt.subplots(5,1,figsize = (7,15))
ax[0].scatter(np.arange(len(lickings1)),lickings1,color = 'k',alpha = 0.4)
count = 0
for i in range(5):
    ax[0].axvline(x = count+trial_nums1[i])
    count += trial_nums1[i]
ax[1].scatter(np.arange(len(lickings2)),lickings2,color = 'k',alpha = 0.4)
count = 0
for i in range(5):
    ax[1].axvline(x = count+trial_nums2[i])
    count += trial_nums2[i]
ax[2].scatter(np.arange(len(lickings3)),lickings3,color = 'k',alpha = 0.4)
count = 0
for i in range(5):
    ax[2].axvline(x = count+trial_nums3[i])
    count += trial_nums3[i]
ax[3].scatter(np.arange(len(lickings4)),lickings4,color = 'k',alpha = 0.4)
count = 0
for i in range(5):
    ax[3].axvline(x = count+trial_nums4[i])
    count += trial_nums4[i]
ax[4].scatter(np.arange(len(lickings5)),lickings5,color = 'k',alpha = 0.4)
count = 0
for i in range(5):
    ax[4].axvline(x = count+trial_nums5[i])
    count += trial_nums5[i]

plt.ylabel('licking rate (/s)')
savepath = 'D:/PhD/Photometry/results/plots/lickings'
plt.savefig("{0}/background_trials_antilicking_during_degradation_{1}.png".format(savepath,date.today()), 
            bbox_inches="tight", dpi = 200)
plt.show()



























