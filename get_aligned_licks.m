function [data_toplot_smoothed,CSin,CSout,TraceStart,TraceEnd,ITIstart] = get_aligned_licks (SR,ThisType,Data,posTrials,time,alignTo,protocol)

%  This function is to be called by the 'plot_cellactivity' function. The
%  input must be the the behavioral+imaging data AFTER it has been reorganized by
%  'DistributionalRL_6Odours_reorganized'. Otherwise it might not work
%  because after this function the data structure is reorganized in a
%  different way from the original saved by Bpod.

plot_auxfigs = 0;

switch alignTo
    case 'CS'
        alignstate = 'Foreperiod';
    case 'US'
        alignstate = 'Trace';
end

if strcmp(protocol,'DRL6Od') || strcmp(protocol,'DRL3Od')
    StimState = ['Stimulus' num2str(ThisType) 'Delivery'];
elseif strcmp(protocol,'VarMgR') || strcmp(protocol,'Eshl16') || strcmp(protocol,'VarMag')
    if ThisType == 2 % this is trial type 4 in this protocol when all trial types are presented, will need to modify to make it more general
        ThisTypeCorr = 7;
    elseif ThisType == 3
        ThisTypeCorr = 2;  % any random StimulusXDelivery that will not ever be called in this protocol
    else ThisTypeCorr = 1;
    end
    StimState = ['Stimulus' num2str(ThisTypeCorr) 'Delivery'];
elseif strcmp(protocol,'Es163O')
    StimState = ['Stimulus' num2str(ThisType) 'Delivery'];
elseif strcmp(protocol,'Revers') || strcmp(protocol,'DRL6Rv')
    StimState = ['Stimulus' num2str(ThisType) 'Delivery'];
elseif strcmp(protocol,'WDMSFP') 
    StimState = 'none';
end


% Parameters for smoothing the lick trace
smooth_win = gausswin(1000,1000/(2*50))/sum(gausswin(1000,1000/(2*50)));  % 100 samples smoothing % 1000/(2*100)= 5
smooth_mult = 1000;

if SR == 1000         % This is the SR I will use in general
    decround = 3;     
elseif SR == 10000 
    decround = 4;     % This is the precision with which the data was acquired (0.1 ms)
elseif SR == 60
    decround = 3;
end    

data_toplot = nan(length(posTrials), length(time));
CSin = nan(1,length(posTrials)); CSout = CSin; ITIstart = CSin;
TraceStart = CSin; TraceEnd = CSin;

for i = 1:length(posTrials)   % Loop across trials of this type  
    
    aligning_time = Data.States(posTrials(i)).(alignstate)(2) + 0.001; % align to the sample 1ms after the end of the foreperiod/trace period (depending on CS or US alignment)        
    trialstart = 0 - aligning_time;
    trialstart_pos = find(round(time,decround) == round(trialstart,decround));
%     trialstart_pos = find(round(time,decround) > round(trialstart,decround),1,'first');
    trialend = Data.States(posTrials(i)).ITI(2) - aligning_time;
    trialend_pos = find(round(time,decround) == round(trialend,decround));
%     trialend_pos = find(round(time,decround) > round(trialend,decround),1,'first');
    if isempty(trialstart_pos); trialstart_pos = 1; end
    if isempty(trialend_pos); trialend_pos = size(data_toplot,2); end
    
    data_toplot(i,trialstart_pos:trialend_pos) = 0;
 
    if isfield(Data.States(posTrials(i)),StimState)
        CSin(i) = Data.States(posTrials(i)).(alignstate)(2) - aligning_time + 0.001;
        CSout(i) = Data.States(posTrials(i)).(StimState)(1,2) - aligning_time + 0.001;
    else                                                                             % This is in case there are trials without CS
        CSin(i) = Data.States(posTrials(i)).Foreperiod(1,2) - aligning_time + 0.001;  
        CSout(i) = CSin(i);
    end
    if isfield(Data.States(posTrials(i)),'Trace')
        TraceStart(i) = Data.States(posTrials(i)).Trace(1,1) - aligning_time + 0.001;
        TraceEnd(i) = Data.States(posTrials(i)).Trace(1,2) - aligning_time + 0.001;
    end
    ITIstart(i) = Data.States(posTrials(i)).ITI(1,1) - aligning_time + 0.001;
    
    if isfield(Data.Events(posTrials(i)),'Port1In')     % Look for the event of interest in this period
        lickin = Data.Events(posTrials(i)).Port1In - aligning_time;
        lickin_pos = nan(size(lickin));
        for j = 1:length(lickin)
            if lickin(j) < time(end)
                lickin_pos(j) = find(round(time,decround) == round(lickin(j),decround));
                data_toplot(i,lickin_pos(j)) = 1;
            end
        end
    end
end

data_toplot_smoothed = smooth_mult*conv2(data_toplot', smooth_win, 'same');

% for k = 1:size(data_toplot,1)
%     if plot_auxfigs
%         figure; hold on
%         plot(time,data_toplot(k,:))
%         plot(time,data_toplot_smoothed(:,k))
%         pause
%     end
% end

if plot_auxfigs
    figure; pcolor(time,1:size(data_toplot,1),data_toplot); shading flat
    hold on; 
    for i = 1:size(data_toplot,1)
        line([CSin(i) CSin(i)],[i i+1],'color',[0 0 0])
        line([CSout(i) CSout(i)],[i i+1],'color',[0 0 0])
        line([TraceEnd(i)+0.001 TraceEnd(i)+0.001],[i i+1],'color',[0 0 0])
        line([ITIstart(i) ITIstart(i)],[i i+1],'color',[0 0 0])
    end
    figure; pcolor(time,1:size(data_toplot_smoothed,1),data_toplot_smoothed); shading flat
    hold on; 
    for i = 1:size(data_toplot,1)
        line([CSin(i) CSin(i)],[i i+1],'color',[0 0 0])
        line([CSout(i) CSout(i)],[i i+1],'color',[0 0 0])
        line([TraceEnd(i)+0.001 TraceEnd(i)+0.001],[i i+1],'color',[0 0 0])
        line([ITIstart(i) ITIstart(i)],[i i+1],'color',[0 0 0])
    end
end
end