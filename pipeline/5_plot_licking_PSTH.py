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
from load_all import Mouse_data
from load_all import pickle_dict
from load_all import load_pickleddata
from sklearn.linear_model import LinearRegression
import scipy.stats as stats
from scipy import signal
from statsmodels.gam.api import GLMGam
import statsmodels.gam.smooth_basis as sb
import statsmodels.api as sm

#%%

mice = ['FgDA_01','FrgD1_01','FrgD2_01','FrgD2_02']
# mice = ['FrgD1_01']
path = 'D:/PhD/Photometry/DATA/round_20220307'

trialtypes = ['go', 'no_go', 'go_omit', 'background','UnpredReward']

licking_data = {}


# get trials

sample_points = 45
for mouse_id in mice:
    load_path = os.path.join(path,'processed/{0}_stats.pickle'.format(mouse_id))
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
save_path = path+'/processed'
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

gap = 10

x = np.linspace(0,1,10)
knots = sb.get_knots_bsplines(x,df = 7)#spacing = 'equal')

basis = sb._eval_bspline_basis(x,degree = 3,knots=knots)[0][:,1]

mice = ['FgDA_01','FrgD1_01','FrgD2_01','FrgD2_02']
# mice = ['FgDA_01']
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
    
        
        savepath = 'D:/PhD/Photometry/DATA/round_20220307/figures'
        plt.savefig("{0}/{1}_{2}_licking.png".format(savepath,mouse_name,TT), bbox_inches="tight", dpi = 72)
        plt.show()


