# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 15:29:38 2022

@author: qianl
"""
import scipy
import numpy as np
import pandas as pd
import matplotlib as plt
from os.path import dirname, join as pjoin
import scipy.io as sio
import math
from datetime import datetime
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



class Mouse_data:
    def __init__(self,mouse_id,protocol,filedir,group = 'T',implant_side = 'L'):
        self.mouse_id = mouse_id
        self.filedir = filedir
        self.filename = ''
        self.protocol = protocol
        self.selected_filename = ''
        self.all_days = []
        self.group = group
        self.training_type = []
        self.df_trials = {}
        self.trialtypes = []
        self.df_trials_iscorrect = {}
        self.df_trials_lick = {}
        self.df_eventcode = {}
        self.p_hit = {}
        self.p_correj = {}
        self.licking_actionwindow = {}
        self.licking_latency = {}
        self.licking_baselicking = {}
        self.stats = {}
        self.event_data = ''
        self.odor_bef = 2.0
        self.odor_on = 1.0
        self.delay = 2.5
        self.rew_after = 8
        self.avg_ITI = 2
        self.implant_side = implant_side
        

    def read_filename(self):
        filedir = pjoin(self.filedir, '{}/{}/Session Data/processed/'.format(self.mouse_id,self.protocol))
    
        filename = []
        for dirpath, dirnames, files in os.walk(filedir): # can walk through all levels down
        #     print(f'Found directory: {dirpath}')
            for f_name in files:
                if f_name.endswith('.mat'):
                    filename.append(dirpath+'/'+f_name)
                    print(f_name)
        print('---------------------------------------------')    
        print('The files have been loaded from the following paths')
        
        self.filename = filename
        

        
    def create_dataset(self): #{'date':df of eventcode} format from original CSV
        date_list = []
        df = {}
        df_2 = {}
        task = []
        
        for session, file in enumerate(self.filename):
            
            perdate_df = {}
            date = re.search(r"(\d{8})",file[-50:-1]).group(0) # extract date: must be like format 2020-02-10
            
            date_list.append(date) # create a list of emperiment date
            
            train_type = os.path.split(file)[-1][len(self.mouse_id)+len(self.protocol)+11:-4]
            
            perdate_df['task'] = train_type
            task.append(train_type) ###
            print(file)
                
            
            mat_contents = sio.loadmat(file) 
            perdate_df['mat_contents'] = mat_contents
            perdate_df['num_ROI'] = mat_contents['MSFPAcq']['ROI'][0,0][0,0]
            perdate_df['num_CAM'] = mat_contents['MSFPAcq']['CAM'][0,0][0,0]
            perdate_df['num_EXC'] = mat_contents['MSFPAcq']['EXC'][0,0][0,0]
            perdate_df['num_FrameRate'] = mat_contents['MSFPAcq']['FrameRate'][0,0][0,0]
            trial_num = mat_contents['Events'].shape[1]
            perdate_df['trial_num'] = trial_num
            perdate_df['session'] = session
            # suppose we keep ROI 0-7,10-11, each list means Time, ROI0, ROI1,...etc
            list_ROIs_gcamp = [[] for _ in range(perdate_df['num_ROI'] +1)]
            list_ROIs_isos = [[] for _ in range(perdate_df['num_ROI'] +1)]
            
            #if perdate_df['num_CAM']== 1 and perdate_df['num_EXC'] == 2:
            for trialnum in range(trial_num):
                for field in range(perdate_df['num_ROI']+1):
                    if field == 0 :
                        list_ROIs_isos[field].append(mat_contents['MSFP']['Time'][0][trialnum][0][0][0]-mat_contents['MSFP']['Time'][0][trialnum][0][0][0][0])
                        list_ROIs_gcamp[field].append(mat_contents['MSFP']['Time'][0][trialnum][0][1][0]-mat_contents['MSFP']['Time'][0][trialnum][0][1][0][0])
                    else: 
                        list_ROIs_isos[field].append(mat_contents['MSFP']['Fluo'][0][trialnum][0][0][:,field-1])
                        list_ROIs_gcamp[field].append(mat_contents['MSFP']['Fluo'][0][trialnum][0][1][:,field-1])
            
    
    
            # ROI correspondence table
            if perdate_df['num_ROI'] == 6 and self.implant_side == 'L':
                corr_ROI = {'AMOT':0,'PMOT':4,'ALOT':1,'PLOT':5,'MNacS':2,'LNacS':3}
                perdate_df['corr_ROI'] = corr_ROI
            elif perdate_df['num_ROI'] == 6 and self.implant_side == 'R':
                corr_ROI = {'AMOT':1,'PMOT':0,'ALOT':4,'PLOT':5,'MNacS':2,'LNacS':3}
                perdate_df['corr_ROI'] = corr_ROI
                
            elif perdate_df['num_ROI'] == 12:
                corr_ROI = {'AMOT':9,'PMOT':8,'ALOT':1,'PLOT':4,'MNacS':5,'LNacS':0,'NacC':2, 'AmLOT':3, 'Pir':10, 'VP':11}
                perdate_df['corr_ROI'] = corr_ROI
            
            if  perdate_df['num_ROI'] == 6:
                
                d = {
                     'TrialStart':mat_contents['States']['TrialStart'][0],
                         'Foreperiod':mat_contents['States']['Foreperiod'][0],
                    'go':mat_contents['States']['go'][0],
                    'no_go':mat_contents['States']['no_go'][0],
                        'go_omit':mat_contents['States']['go_omit'][0],
                        'c_odor':mat_contents['States']['c_odor'][0],
                        'c_odor_omit':mat_contents['States']['c_odor_omit'][0],
                    'background':mat_contents['States']['background'][0],
                    'Trace':mat_contents['States']['Trace'][0],
                    'ITI':mat_contents['States']['ITI'][0],
                    'UnpredReward':mat_contents['States']['UnexpectedReward'][0],
                    'water':mat_contents['States']['Reward'][0],
                    'TrialEnd':mat_contents['States']['TrialEnd'][0],
                     'AMOT': list_ROIs_gcamp[corr_ROI['AMOT']+1], 'PMOT': list_ROIs_gcamp[corr_ROI['PMOT']+1], 'ALOT':list_ROIs_gcamp[corr_ROI['ALOT']+1],
                             'PLOT':list_ROIs_gcamp[corr_ROI['PLOT']+1],'MNacS':list_ROIs_gcamp[corr_ROI['MNacS']+1],'LNacS':list_ROIs_gcamp[corr_ROI['LNacS']+1],
                             'Doric_Time_EXC1':list_ROIs_gcamp[0],
                             'AMOT_isos': list_ROIs_isos[corr_ROI['AMOT']+1], 'PMOT_isos': list_ROIs_isos[corr_ROI['PMOT']+1], 
                             'ALOT_isos':list_ROIs_isos[corr_ROI['ALOT']+1],'PLOT_isos':list_ROIs_isos[corr_ROI['PLOT']+1],
                             'MNacS_isos':list_ROIs_isos[corr_ROI['MNacS']+1],
                             'LNacS_isos':list_ROIs_isos[corr_ROI['LNacS']+1],
                             }
                
            elif perdate_df['num_ROI'] == 12:
                d = {
                     'TrialStart':mat_contents['States']['TrialStart'][0],
                         'Foreperiod':mat_contents['States']['Foreperiod'][0],
                    'go':mat_contents['States']['go'][0],
                    'no_go':mat_contents['States']['no_go'][0],
                        'go_omit':mat_contents['States']['go_omit'][0],
                        'c_odor':mat_contents['States']['c_odor'][0],
                        'c_odor_omit':mat_contents['States']['c_odor_omit'][0],
                    'background':mat_contents['States']['background'][0],
                    'Trace':mat_contents['States']['Trace'][0],
                    'ITI':mat_contents['States']['ITI'][0],
                    'UnpredReward':mat_contents['States']['UnexpectedReward'][0],
                    'water':mat_contents['States']['Reward'][0],
                    'TrialEnd':mat_contents['States']['TrialEnd'][0],
                     'AMOT': list_ROIs_gcamp[corr_ROI['AMOT']+1], 'PMOT': list_ROIs_gcamp[corr_ROI['PMOT']+1], 'ALOT':list_ROIs_gcamp[corr_ROI['ALOT']+1],
                             'PLOT':list_ROIs_gcamp[corr_ROI['PLOT']+1],'MNacS':list_ROIs_gcamp[corr_ROI['MNacS']+1],'LNacS':list_ROIs_gcamp[corr_ROI['LNacS']+1],
                             'NacC':list_ROIs_gcamp[corr_ROI['NacC']+1],'AmLOT':list_ROIs_gcamp[corr_ROI['AmLOT']+1],
                             'Pir':list_ROIs_gcamp[corr_ROI['Pir']+1],'VP':list_ROIs_gcamp[corr_ROI['VP']+1],
                             
                             'Doric_Time_EXC1':list_ROIs_gcamp[0],
                             'AMOT_isos': list_ROIs_isos[corr_ROI['AMOT']+1], 'PMOT_isos': list_ROIs_isos[corr_ROI['PMOT']+1], 
                             'ALOT_isos':list_ROIs_isos[corr_ROI['ALOT']+1],'PLOT_isos':list_ROIs_isos[corr_ROI['PLOT']+1],
                             'MNacS_isos':list_ROIs_isos[corr_ROI['MNacS']+1],
                             'LNacS_isos':list_ROIs_isos[corr_ROI['LNacS']+1],
                             'NacC_isos':list_ROIs_isos[corr_ROI['NacC']+1],'AmLOT_isos':list_ROIs_isos[corr_ROI['AmLOT']+1],
                             'Pir_isos':list_ROIs_isos[corr_ROI['Pir']+1],'VP_isos':list_ROIs_isos[corr_ROI['VP']+1],
                             }
            d_2 = {
                 'TrialStart':[x[0] for x in mat_contents['States']['TrialStart'][0]],
                     'Foreperiod':[x[0] for x in mat_contents['States']['Foreperiod'][0]],
                'go_odor':[x[0] for x in mat_contents['States']['go'][0]],
                'nogo_odor':[x[0] for x in mat_contents['States']['no_go'][0]],
                    'go_omit':[x[0] for x in mat_contents['States']['go_omit'][0]],
                    'control_odor':[x[0] for x in mat_contents['States']['c_odor'][0]],
                    'control_odor_omit':[x[0] for x in mat_contents['States']['c_odor_omit'][0]],
                'background':mat_contents['States']['background'][0],

                'UnpredReward':[x[0] for x in mat_contents['States']['UnexpectedReward'][0]],
                'water_on':[x[0][0] for x in mat_contents['States']['Reward'][0]],
                'water_off':[x[0][1] for x in mat_contents['States']['Reward'][0]],
                'trial_end':[x[0][1] for x in mat_contents['States']['TrialEnd'][0]]}
                 
            
            df_trial = pd.DataFrame(data = d)  
            df_trial_2 = pd.DataFrame(data = d_2)  
            
            # add trialtype to the dataframe
            trialtype = []
            for index, row in df_trial_2.iterrows():
                if not math.isnan(row['go_odor'][0]):
                    trialtype.append('go')
                elif not math.isnan(row['nogo_odor'][0]):
                    trialtype.append('no_go')
                elif not math.isnan(row['go_omit'][0]):
                    trialtype.append('go_omit')
                    df_trial_2.at[index, 'go_odor'] =row['go_omit']
                elif not math.isnan(row['background'][0][0]):
                    trialtype.append('background')
                elif not math.isnan(row['UnpredReward'][0]):
                    trialtype.append('unpred_water')
                elif not math.isnan(row['control_odor'][0]):
                    trialtype.append('c_reward')
                elif not math.isnan(row['control_odor_omit'][0]):
                    trialtype.append('c_omit')
                    df_trial_2.at[index, 'control_odor'] =row['control_odor_omit']

            
            df_trial.insert(0,'Trialtype',trialtype)
            df_trial_2.insert(0,'trialtype',trialtype)
            df_trial['lickings'] = [x[0] if len(x) != 0 else x for x in mat_contents['Events']['Port1Out'][0]]
            df_trial_2['licking'] = [x[0] if len(x) != 0 else x for x in mat_contents['Events']['Port1Out'][0]]
            df_trial['id'] = [self.mouse_id]*trial_num
            df_trial['trialnum'] = np.arange(trial_num)
            df_trial['session'] = [session] *trial_num
            df_trial['task'] = [perdate_df['task']] *trial_num
            df_trial['group'] = [self.group] *trial_num
            
            perdate_df['dataframe'] = df_trial

            
            df.update({str(session) +'_'+ date:perdate_df}) # create the dict of key: date and value: data dataframe
            
            df_2.update({str(session) +'_'+ date:df_trial_2})
        self.df_bpod_doric = df #individual mouse event code data
        self.df_trials = df_2
        date_format = '%Y-%m-%d'
        index = np.argsort(date_list)
        self.all_days = [str(session) + '_' + date_list[i] for i in index]
        self.training_type = [task[i] for i in index]
        
        print('---------------------------------------------')
        print('{0} has data from these days: {1}'.format(self.mouse_id,list(zip(self.all_days,self.training_type))))
        


def pickle_dict(df,path,filename):
    try:
        os.makedirs(path) # create the path first
    except FileExistsError:
        print('the path exist.')
    filename = path +'/{}.pickle'.format(filename)
    with open(filename, 'wb') as handle:
        pickle.dump(df, handle, protocol=pickle.HIGHEST_PROTOCOL)
    print('save to pickle done!')


def load_pickleddata(filename):
    
    with open(filename, 'rb') as handle:
        df = pickle.load(handle)
    return df

#%% main code
if __name__ == '__main__':
    
    is_save = True

    
    #********************
    load_path = 'D:/PhD/Photometry/DATA/round_20221111-right-pilot_c-odor'
    
    # load file
  
    
    # mouse_names = ['FgDA_03','FgDA_04','FgDA_05']
    # mouse_names = ['FgDA_02']
    mouse_names = ['FgDA_07']
    for mouse_name in mouse_names:
        cute = Mouse_data(mouse_name, protocol = 'Selina_C5D5R3E5R3',filedir = load_path,group = 'C',implant_side = 'R') #### group 
        cute.read_filename()
        #parse data
        cute.create_dataset()
        
        if mouse_name == 'FgDA_01':
            cute.df_bpod_doric['3_20220214']['dataframe'] = cute.df_bpod_doric['3_20220214']['dataframe'].iloc[0:80]
        elif mouse_name == 'FgDA_02':

            cute.df_bpod_doric['6_20220511']['dataframe'] = cute.df_bpod_doric['6_20220511']['dataframe'].iloc[0:50]


        elif mouse_name == 'FgDA_03':
            cute.df_bpod_doric['5_20220509']['dataframe'] = cute.df_bpod_doric['5_20220509']['dataframe'].iloc[0:70]
            cute.df_bpod_doric['6_20220510']['dataframe'] = cute.df_bpod_doric['6_20220510']['dataframe'].iloc[10:110]
            cute.df_bpod_doric['7_20220511']['dataframe'] = cute.df_bpod_doric['7_20220511']['dataframe'].iloc[20:120]
            cute.df_bpod_doric['8_20220512']['dataframe'] = cute.df_bpod_doric['8_20220512']['dataframe'].iloc[0:60]
            cute.df_bpod_doric['9_20220513']['dataframe'] = cute.df_bpod_doric['9_20220513']['dataframe'].iloc[0:60]
        elif mouse_name == 'FgDA_05':
            cute.df_bpod_doric['1_20220503']['dataframe'] = cute.df_bpod_doric['1_20220503']['dataframe'].iloc[20:160]
            cute.df_bpod_doric['2_20220504']['dataframe'] = cute.df_bpod_doric['2_20220504']['dataframe'].iloc[0:80]
            cute.df_bpod_doric['3_20220506']['dataframe'] = cute.df_bpod_doric['3_20220506']['dataframe'].iloc[0:60]
            cute.df_bpod_doric['4_20220508']['dataframe'] = cute.df_bpod_doric['4_20220508']['dataframe'].iloc[20:120]
            cute.df_bpod_doric['5_20220509']['dataframe'] = cute.df_bpod_doric['5_20220509']['dataframe'].iloc[20:80]
            cute.df_bpod_doric['6_20220510']['dataframe'] = cute.df_bpod_doric['6_20220510']['dataframe'].iloc[0:80]
            cute.df_bpod_doric['10_20220516']['dataframe'] = cute.df_bpod_doric['10_20220516']['dataframe'].iloc[0:140]
            cute.df_bpod_doric['12_20220518']['dataframe'] = cute.df_bpod_doric['12_20220518']['dataframe'].iloc[40:140]
            cute.df_bpod_doric['13_20220519']['dataframe'] = cute.df_bpod_doric['13_20220519']['dataframe'].iloc[0:40]
            cute.df_bpod_doric['18_20220525']['dataframe'] = cute.df_bpod_doric['18_20220525']['dataframe'].iloc[20:120]
                           
            
            
        if is_save:
            #save data by pickle
            #****************
            save_path = load_path+'/processed'
            
            filename = '{}_stats'.format(cute.mouse_id)
            pickle_dict(cute,save_path,filename)
            
   