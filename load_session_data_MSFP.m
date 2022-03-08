% function [licks_alltrials, time_licks, fluo_output_alltrials, time_fluo, ...
%     revtrial, TT_allRsizes, thisRsize_pos, CSin, CSout,...
%     TraceEnd,SR_fluo,num_EXC, RotEnc_output_alltrials] = load_session_data_MSFP(protocol, mouse, day, alignTo,...
%     sinal_filtLicking)

% Load the data
% pathhere = mfilename('fullpath');
% preproc_data_folder = [pathhere(1:end-22) '/AllData'];
% preproc_data_folder = '/home/sara/Documents/DATA/MultiSiteFiberPhotometry/Analysis/AllData';
% 
% load(fullfile(preproc_data_folder, [protocol '_' mouse '_' day '.mat']),'Session'); 
% load(fullfile(preproc_data_folder, [protocol '_' mouse '_' day '.mat']),'States');  
% load(fullfile(preproc_data_folder, [protocol '_' mouse '_' day '.mat']),'Events');   
% % load(fullfile(preproc_data_folder, [protocol '_' mouse '_' day '.mat']),'RotEncoder'); 
% load(fullfile(preproc_data_folder, [protocol '_' mouse '_' day '.mat']),'MSFP');
% load(fullfile(preproc_data_folder, [protocol '_' mouse '_' day '.mat']),'MSFPAcq');

Data.States = States;
Data.Events = Events;
Data.Session = Session;
% Data.RotEncoder = RotEncoder;
Data.MSFP = MSFP;

% Sampling rates to use:
SR_licks = 1000;
time_licks = -6:1/SR_licks:20-1/SR_licks;
SR_fluo = mean(MSFPAcq.FrameRate);
time_fluo = -6:1/SR_fluo:20-1/SR_fluo;
SR_RotEnc = [];  % the rotary encoder does not have a fixed sampling rate: it is acquired when movement is detected
time_RotEnc = time_fluo; % I will have these two in the same time to facilitate analysis (correlations, etc), so the Rotary Encoder will be downsampled

% Trial types
TT = unique(Session.TrialTypes);   
num_TT = length(TT);
num_trials = length(MSFP);

% Fiber Photometry acquisition
num_EXC = MSFPAcq.EXC;
num_rois = MSFPAcq.ROI;

sinal_filtLicking = 1;
alignTo = 'CS';
protocol = 'Selina_C5D5R3E3R3';
%% See if we should ignore trials after mouse looses motivation

if sinal_filtLicking
      
    posAllTrials = 1:num_trials;
    
    % Get aligned licks
    [licks_all,CSin,~,~, TraceEnd,...
        ~] =  get_aligned_licks(SR_licks,nan,...
        Data,posAllTrials,time_licks,alignTo, protocol);
    
    % mean lick to CS
    CSstart_samples = find(time_licks>mode(CSin),1,'first');               % CS quantification start
    TraceEnd_samples = find(time_licks<mode(TraceEnd),1,'last');           % Trace quantification start
    meanlick = nanmean(licks_all(CSstart_samples:TraceEnd_samples,:),1);
    
    limit_trials = find(movmean(meanlick(21:end),15)==0,1,'first');
    if isempty(limit_trials); limit_trials = nan; end
    limit_trials = limit_trials + 20;
    
    figure; hold on
    plot(meanlick)
    plot(movmean(meanlick,10))
    line([limit_trials limit_trials],[0 10],'color','r')
    title([mouse '-' day]) 
else
    limit_trials = num_trials;
end

%% Get all fluo info for all trials, separated by trial type.  
posTrials = cell(1,num_TT);
fluo_output_alltrials = cell(1,num_TT);
licks_alltrials = cell(1,num_TT);
CSin_all = cell(1,num_TT);
CSout_all = cell(1,num_TT);
TraceStart_all = cell(1,num_TT);
TraceEnd_all = cell(1,num_TT);
ITIstart_all = cell(1,num_TT);
revtrial = nan(1,num_TT);
TT_allRsizes = cell(1,num_TT);

for i_TT = 1:num_TT
   posTrials{i_TT} = find(Session.TrialTypes == TT(i_TT));
   posTrials{i_TT} = posTrials{i_TT}(posTrials{i_TT} <= nanmin(num_trials,limit_trials));  % To ignore trials that were not imaged (due to some issue) (rare condition) or that happened after mouse lost motivation
   
   [licks_alltrials{i_TT},CSin_all{i_TT},CSout_all{i_TT},...
       TraceStart_all{i_TT}, TraceEnd_all{i_TT},...
       ITIstart_all{i_TT}] = get_aligned_licks(SR_licks,i_TT,...
       Data,posTrials{i_TT},time_licks,alignTo, protocol); 
   % licks_alltrials is a matrix with time in rows (time_licks at with 1ms precision) and trials in columns
   
   fluo_output_alltrials{i_TT} = get_aligned_fluo_MSFP(num_rois,...
       TT(i_TT),Data,posTrials{i_TT},time_fluo,alignTo);
   % dimensions of fluo_output: time,trials,rois
   
   RotEnc_output_alltrials{i_TT} = get_aligned_RotEnc(TT(i_TT),...
       Data,posTrials{i_TT},time_RotEnc,alignTo,protocol);
   
   % This will only be used in the Reversal Task
   if strcmp(protocol,'Revers') 
       outcomes_thistriatltype = Session.RewardDelivered(posTrials{i_TT});
       if ~isempty(find(diff(outcomes_thistriatltype)~=0,1,'first'))
           revtrial(i_TT) = find(diff(outcomes_thistriatltype)~=0,1,'first')+1;  % First trial with reversed outcome for this CS
       else
           revtrial(i_TT) = nan;
       end
   elseif strcmp(protocol,'DRL6Rv')
       outcomes_thistriatltype = Session.RewardDelivered(posTrials{i_TT});
       if i_TT < 4 && ~isempty(find(outcomes_thistriatltype ~= 0,1))
           revtrial(i_TT) = find(outcomes_thistriatltype ~= 0,1,'first');  % First trial with reversed outcome in this protocol is when reward is larger than zero for TT1-3
       elseif i_TT >= 4 && ~isempty(find(outcomes_thistriatltype == 0,1))
           revtrial(i_TT) = find(outcomes_thistriatltype == 0,1,'first');  % First trial with reversed outcome in this protocol is when reward is larger than zero for TT4-6
       else
           revtrial(i_TT) = nan;
       end
   else      
       revtrial(i_TT) = nan;
   end
   
   % Get reward sizes and location for each trial type
   if isfield(Session,'RewardDelivered')
        RewDelivName = 'RewardDelivered';
   elseif isfield(Session,'Rewards')
        RewDelivName = 'Rewards';
   end
   TT_allRsizes{i_TT} = unique(Session.(RewDelivName)(posTrials{i_TT}));
   for n = 1:length(TT_allRsizes{i_TT})
       thisRsize_pos{i_TT,n} = find(Session.(RewDelivName)(posTrials{i_TT}) == TT_allRsizes{i_TT}(n));
   end
end

% In the data I am analysing it will the the same for all trials in a session
CSin = CSin_all{1}(1);
CSout = CSout_all{1}(1);
TraceEnd = TraceEnd_all{1}(1);

% end