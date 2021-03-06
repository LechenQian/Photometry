function [output] = get_aligned_fluo_MSFP(num_rois,ThisType,TrialData,posTrials,time,alignTo)



dbstop if error
% The DF0F is done based on F0 being just before the CS start. For
% quantification of the US, I will later subtract tha mean value a few
% 100-200ms before US start

plot_auxfigs = 0;

%% Define aligning state and output type

switch alignTo
    case 'CS'
        alignstate = 'Foreperiod';
    case 'US'
        alignstate = 'Trace';
end

calcDFoF = 1;           % in this case, calculate DFoF with F0 mean just before CS start
num_EXC = size(TrialData.MSFP(1).Time,2);

%% Align and calculate DFoF if necessary

output = nan(length(time),length(posTrials),num_rois,num_EXC);
figure;hold on
for i = 1:length(posTrials)   % Loop across trials of this type
    
    aligning_time = TrialData.States(posTrials(i)).(alignstate)(2) + 0.001; % align to the sample 1ms after the end of the foreperiod/trace period (depending on CS or US alignment)
    ForeperiodDur = TrialData.States(posTrials(i)).(alignstate)(2) - TrialData.States(posTrials(i)).(alignstate)(1)+0.001;
    if num_EXC == 3
        if ~isempty(TrialData.MSFP(posTrials(i)).Time{1}) && ~isempty(TrialData.MSFP(posTrials(i)).Time{2}) && ~isempty(TrialData.MSFP(posTrials(i)).Time{3})
            
            for i_exc = 1:num_EXC
    
                exctrialtime = TrialData.MSFP(posTrials(i)).Time{i_exc}-TrialData.MSFP(posTrials(i)).Time{i_exc}(1);
                aligned_exctrialtime = exctrialtime - ForeperiodDur; % add the time of the Foreperiod start, because the TTL goes up at the start of Foreperiod.
                
                % F0 window
                F0_start = find(aligned_exctrialtime > -1.0, 1, 'first'); % Originally I had it startin at -1.0, for analysis of baseline flutuations I might want to use from -0.6
                F0_end = find(aligned_exctrialtime < -0.2, 1, 'last');
    
                % Fluo trace from different origins
                fluo_traces = TrialData.MSFP(posTrials(i)).Fluo{i_exc};  % rows are fluo values in time, columns are ROIs
    
                if calcDFoF
                    baseline = repmat(mean(fluo_traces(F0_start:F0_end,:)),length(aligned_exctrialtime),1);
                    fluo_traces_outuput = (fluo_traces - baseline)./baseline;
                else
                    fluo_traces_outuput = fluo_traces - mean(fluo_traces(F0_start:F0_end,:)); 
                end    
                
                if ThisType ~= 1 && i_exc == 2
                    plot(aligned_exctrialtime,fluo_traces_outuput(:,5)+10*i)
                end
                % Resample to have the same samples size across trials:
                % Samples used to delimit the resampling
                pos1 = find(time-aligned_exctrialtime(1) >= 0, 1);
                pos2 = find(time-aligned_exctrialtime(end) <= 0, 1,'last'); 
    
                serietempo = timeseries(fluo_traces_outuput,aligned_exctrialtime,'Name','Bckgnd');
                res_serietempo = resample(serietempo,time(pos1:pos2));
                output(pos1:pos2,i,:,i_exc) = squeeze(res_serietempo.data); 
            end
        end
    elseif num_EXC == 2
        if ~isempty(TrialData.MSFP(posTrials(i)).Time{1}) && ~isempty(TrialData.MSFP(posTrials(i)).Time{2})
            
            for i_exc = 1:num_EXC
    
                exctrialtime = TrialData.MSFP(posTrials(i)).Time{i_exc}-TrialData.MSFP(posTrials(i)).Time{i_exc}(1);
                aligned_exctrialtime = exctrialtime - ForeperiodDur; % add the time of the Foreperiod start, because the TTL goes up at the start of Foreperiod.
                
                % F0 window
                F0_start = find(aligned_exctrialtime > -1.0, 1, 'first'); % Originally I had it startin at -1.0, for analysis of baseline flutuations I might want to use from -0.6
                F0_end = find(aligned_exctrialtime < -0.2, 1, 'last');
    
                % Fluo trace from different origins
                fluo_traces = TrialData.MSFP(posTrials(i)).Fluo{i_exc};  % rows are fluo values in time, columns are ROIs
    
                if calcDFoF
                    baseline = repmat(mean(fluo_traces(F0_start:F0_end,:)),length(aligned_exctrialtime),1);
                    fluo_traces_outuput = (fluo_traces - baseline)./baseline;
                else
                    fluo_traces_outuput = fluo_traces - mean(fluo_traces(F0_start:F0_end,:)); 
                end    
                
                if ThisType ~= 1 && i_exc == 2
                    plot(aligned_exctrialtime,fluo_traces_outuput(:,5)+10*i)
                end
                % Resample to have the same samples size across trials:
                % Samples used to delimit the resampling
                pos1 = find(time-aligned_exctrialtime(1) >= 0, 1);
                pos2 = find(time-aligned_exctrialtime(end) <= 0, 1,'last'); 
    
                serietempo = timeseries(fluo_traces_outuput,aligned_exctrialtime,'Name','Bckgnd');
                res_serietempo = resample(serietempo,time(pos1:pos2));
                output(pos1:pos2,i,:,i_exc) = squeeze(res_serietempo.data); 
            end
        end
    end
end  


if plot_auxfigs 
    for j = 1 %1:num_rois
        figure; hold on
        subplot(1,2,1); hold on
        pcolor(time,1:size(output,2),squeeze(output(:,:,j,1))'); shading flat
        line([0 0],[1 size(output,2)],'color','r')
        line([3 3],[1 size(output,2)],'color','r')
        line([-1 -1],[1 size(output,2)],'color','r')
        set(gca, 'TickDir', 'out')
        title([' Trials CS' num2str(ThisType) '-ROI ' num2str(j)])
        subplot(1,2,2); hold on
        pcolor(time,1:size(output,2),squeeze(output(:,:,j,2))'); shading flat
        line([0 0],[1 size(output,2)],'color','r')
        line([3 3],[1 size(output,2)],'color','r')
        line([-1 -1],[1 size(output,2)],'color','r')
        set(gca, 'TickDir', 'out')
        title([' Trials CS' num2str(ThisType) '-ROI ' num2str(j)])
        
    end
    
    % To see the fluo trials sequentially:
    figure; hold on; 
    for i=1:length(posTrials)
        plot(TrialData.MSFP(i).Time{2},TrialData.MSFP(i).Fluo{2}(:,5));
    end
    figure; hold on; 
    for i=1:length(posTrials)
        ForeperiodDur = TrialData.States(posTrials(i)).(alignstate)(2) - TrialData.States(posTrials(i)).(alignstate)(1)+0.001;
        plot(TrialData.MSFP(i).Time{2}-TrialData.MSFP(i).Time{2}(1)-ForeperiodDur,TrialData.MSFP(i).Fluo{2}(:,2)+i*100);
    end
    
end
end