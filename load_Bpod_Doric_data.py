import scipy
import numpy as np
import pandas as pd
import matplotlib as plt
from os.path import dirname, join as pjoin
import scipy.io as sio
import math


class Mouse_data:
    def __init__(self,mouse_id,filedir):
        self.mouse_id = mouse_id
        self.filedir = filedir
        self.filename = ''
        self.selected_filename = ''
        self.all_days = []
        
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
        self.odor_bef = 3.0
        self.odor_on = 1.0
        self.delay = 2.5
        self.rew_after = 8
    # input
    mouse_id = 'Test'
    protocol = 'Selina_C5D5R3E5R3'
    
    data_dir = pjoin('D:/PhD/Photometry/DATA/round_20220307/', '{}/{}/Session Data/processed/'.format(mouse_id,protocol))
    
    mat_fname = pjoin(data_dir, 'FgDA_01_Selina_C5D5R3E5R3_20220207_cond1.mat')
    mat_contents = sio.loadmat(mat_fname)
    task = mat_fname[len(data_dir)+len(mouse_id)+4+len(protocol)+10:-4]
    # get fluo data: 4 index, 1st be 0, 2nd trial num, 3rd be 0, 4th exc num, get a matri of signal for ROIs
    # e.g. mat_contents['MSFP']['Fluo'][0][1][0][0] -> (203,12)
    #%%
    
    num_ROI = mat_contents['MSFPAcq']['ROI'][0,0][0,0]
    num_CAM = mat_contents['MSFPAcq']['CAM'][0,0][0,0]
    num_EXC= mat_contents['MSFPAcq']['EXC'][0,0][0,0]
    num_FrameRate = mat_contents['MSFPAcq']['FrameRate'][0,0][0,0]
    trial_num = mat_contents['Events'].shape[1]
    session = 0
    # suppose we keep ROI 0-7,10-11, each list means Time, ROI0, ROI1,...etc
    list_ROIs_gcamp = [[] for _ in range(num_ROI +1)]
    list_ROIs_isos = [[] for _ in range(num_ROI +1)]
    
    if num_CAM == 1 and num_EXC == 2:
        for trialnum in range(trial_num):
            for field in range(num_ROI+1):
                if field == 0 :
                    list_ROIs_gcamp[field].append(mat_contents['MSFP']['Time'][0][trialnum][0][0][0]-mat_contents['MSFP']['Time'][0][trialnum][0][0][0][0])
                    list_ROIs_isos[field].append(mat_contents['MSFP']['Time'][0][trialnum][0][1][0]-mat_contents['MSFP']['Time'][0][trialnum][0][1][0][0])
                else: 
                    list_ROIs_gcamp[field].append(mat_contents['MSFP']['Fluo'][0][trialnum][0][0][:,field-1])
                    list_ROIs_isos[field].append(mat_contents['MSFP']['Fluo'][0][trialnum][0][1][:,field-1])
    elif num_CAM == 2 and num_EXC == 3:
        pass
    
    
    # ROI correspondence table
    if num_ROI == 6:
        corr_ROI = {'AMOT':0,'PMOT':4,'ALOT':5,'PLOT':3,'MNacS':2,'LNacS':3}
    elif num_ROI == 12:
        corr_ROI = {'AMOT':9,'PMOT':8,'ALOT':1,'PLOT':4,'MNacS':5,'LNacS':0,'NacC':2, 'AmLOT':3, 'Pir':10, 'VP':11}
                
    d = {
         'TrialStart':mat_contents['States']['TrialStart'][0],
             'Foreperiod':mat_contents['States']['Foreperiod'][0],
        'go':mat_contents['States']['go'][0],
        'no_go':mat_contents['States']['no_go'][0],
            'go_omit':mat_contents['States']['go_omit'][0],
        'background':mat_contents['States']['background'][0],
        'Trace':mat_contents['States']['Trace'][0],
        'ITI':mat_contents['States']['ITI'][0],
        'UnpredReward':mat_contents['States']['UnexpectedReward'][0],
        'Reward':mat_contents['States']['Reward'][0],
        'TrialEnd':mat_contents['States']['TrialEnd'][0],
         'AMOT': list_ROIs_gcamp[corr_ROI['AMOT']+1], 'PMOT': list_ROIs_gcamp[corr_ROI['PMOT']+1], 'ALOT':list_ROIs_gcamp[corr_ROI['ALOT']+1],
                 'PLOT':list_ROIs_gcamp[corr_ROI['PLOT']+1],'MNacS':list_ROIs_gcamp[corr_ROI['MNacS']+1],'LNacS':list_ROIs_gcamp[corr_ROI['LNacS']+1],
                 'Doric_Time_EXC1':list_ROIs_gcamp[0],
                 'AMOT_isos': list_ROIs_isos[corr_ROI['AMOT']+1], 'PMOT_isos': list_ROIs_isos[corr_ROI['PMOT']+1], 
                 'ALOT_isos':list_ROIs_isos[corr_ROI['ALOT']+1],'PLOT_isos':list_ROIs_isos[corr_ROI['PLOT']+1],
                 'MNacS_isos':list_ROIs_isos[corr_ROI['MNacS']+1],
                 'LNacS_isos':list_ROIs_isos[corr_ROI['LNacS']+1],
                 }
    df_trial = pd.DataFrame(data = d)  
    
    # add trialtype to the dataframe
    trialtype = []
    for index, row in df_trial.iterrows():
        if not math.isnan(row['go'][0][0]):
            trialtype.append('go')
        elif not math.isnan(row['no_go'][0][0]):
            trialtype.append('no_go')
        elif not math.isnan(row['go_omit'][0][0]):
            trialtype.append('go_omit')
        elif not math.isnan(row['background'][0][0]):
            trialtype.append('background')
        elif not math.isnan(row['UnpredReward'][0][0]):
            trialtype.append('UnpredReward')
    df_trial.insert(0,'Trialtype',trialtype)
    df_trial['lickings'] = mat_contents['Events']['Port1Out'][0]
    df_trial['id'] = [mouse_id]*trial_num
    df_trial['trialnum'] = np.arange(trial_num)
    df_trial['session'] = [session] *trial_num
    df_trial['task'] = [task] *trial_num
    df_trial['group'] = ['T'] *trial_num
    
            
        