function [ind_BSL, ind_CS, ind_US_exp, ind_US_unexp, ind_US_exp_early] = get_indicesEvents(time_fluo,fluotype,CSin,CSstartT,CSendT,TraceEnd)

% CSstartT CS start time: to choose differently for licking and fluo, for
% instance. CSendT same thing for the end of CS quantification period

if strcmp(fluotype,'NPX_ephys')
    USlim1 = 0.2; USlim2 = 0.6;
    
elseif length(fluotype)>15 && strcmp(fluotype(1:16),'Deconv_SpikeRate') % this is based on spike count/bin at 60Hz
    USlim1 = 0.2; USlim2 = 0.6;
    
elseif length(fluotype)>15 && strcmp(fluotype(1:18),'Deconv_spikesTimes') % each spike is convolved with a gamma function
    USlim1 = 0.35; USlim2 = 0.6; %USlim2 = 1.2; this is for convolving with 200 ms gamma
    USlim1 = 0.25; USlim2 = 0.6;
%     USlim1 = 0.20; USlim2 = 0.6; % USlim2 = 0.45 these limits also do not work bad
elseif strcmp(fluotype,'FiberPhotometry')
    USlim1 = 1.0; USlim2 = 1.9;
    
else % For fluo values in 2P acquisitions
    USlim1 = 0.6; USlim2 = 1.5;
end

ind_BSL = find(time_fluo >= CSin{1}-3.0,1,'first'):find(time_fluo <= CSin{1}-0.7,1,'last');
ind_CS = find(time_fluo >= CSin{1}+CSstartT,1,'first'):find(time_fluo <= CSin{1}+CSendT,1,'last');
ind_US_exp_early = find(time_fluo >= TraceEnd{1}+0.01,1,'first'):find(time_fluo <= TraceEnd{1}+0.8,1,'last');
ind_US_exp = find(time_fluo >= TraceEnd{1}+USlim1,1,'first'):find(time_fluo <= TraceEnd{1}+USlim2,1,'last');
ind_US_unexp = find(time_fluo >= CSin{1}+USlim1,1,'first'):find(time_fluo <= CSin{1}+USlim2,1,'last');
end