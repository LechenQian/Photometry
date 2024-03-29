function sync_bpod_doric_data_for_20221006andafter(filedir, filename,DoricStudioVersion)

    
    dbstop if error
    
    % e.g. protocol: 'DistributionalRL_6Odours',
    % 'VariableMagnitudeReward','Eshel2016','ReversalTask',
    % 'DistRL_Rev_6Odors', 'WaterDelivery_MSFP'
    % e.g. mouse: 'SM79'
    % e.g. day: '20190404'
    % e.g. protocolname_save: 'DRL6Od','VarMgR','Eshl16','Revers','DRL6Rv','WDMSFP','Es163O'
    % DoricStudioVersion is the software version used for recording, because in
    % version 5.4.1.23 there is a mismatch in the labels of the headers:CAM1 = CAM1
    % CAM2 = EXC1
    % EXC1 = EXC2
    % EXC2 = EXC3
    % EXC3 = CAM2
    
    % MSFP = multi-site fiber photometry
    
    % S. Matias, November 2021
    % Data will be saved in a structure containing both behavioral and MSFP
    % data called 'Data'. This will later be used to coordinate with video
    % data.
    
    % Quick notes about movement artifact correction
    % See also https://github.com/leomol/FPA
    
    % Options for movement artifact correction:
    % use red channel (tdTomato) or isosbestic wavelenght (415 nm)
    % use the entire trace, or select only ITI period
    
    plotauxfig = 0;
    
    %%  Currently we have no video data (from cameras looking at the mice) - add later
    
    
    %% Load the data saved by the Doric Studio software (from the text file .csv)
    % I will load the .csv file and look at the output of the ROI in the
    % different cameras/acquisition and test the movement artifact correction
    % method that is more appropriate
    
    doric_data_folder = strcat(filedir,'/doric/');
    doric_filename = [strrep(filename,'mat','csv')];
    
    % Get the names of the variables that describe each column (two first lines
    % of the file)
    % Note: the labels might not be always the exact same in all .csv files,
    % depending on the DI/O channel was selected to be included, for instance
    
    cd(doric_data_folder);     % Returns all the files and folders in the directory                  
    doric_file = dir([doric_filename{1}]); 
    
    if ~isempty(doric_file)
    
        fluodataAvailable = 1;
        fprintf('Processing session %s\n',doric_filename{1});
        
        % Load headers on the first row
        fid = fopen([doric_data_folder doric_filename{1}], 'r');
        tmp = fgetl(fid);
        tmp = textscan(tmp, '%s', 'Delimiter', ',');
        tmp = tmp{1};
        names1 = tmp;
        % Load headers on the second row
        tmp = fgetl(fid);
        tmp = textscan(tmp, '%s', 'Delimiter', ',');
        tmp = tmp{1};
        names2 = tmp;
    
        % Combine the headers of the two first lines to create unique variable names   
        num_var = numel(names1);         % total number of variables
        var_names = cell(num_var,1);
        
        
        %% Number of items for analysis
    
        % Number of cameras
        indCAM = strfind(names1, 'CAM');
        IndexCAM = find(not(cellfun('isempty',indCAM)));
        num_CAM = length(IndexCAM);
    
        % Number of excitations
        indEXC = strfind(names1, 'EXC');
        IndexEXC = find(not(cellfun('isempty',indEXC)));
        num_EXC = length(IndexEXC);
    
        % Number of ROIs
        indROI = strfind(names1, 'ROI');
        num_ROI_acquisitions = length(find(not(cellfun('isempty',indROI))));
        num_ROI = num_ROI_acquisitions/num_EXC;
    
        % Save this information in the data strucure
        Data.MSFPAcq.CAM = num_CAM;
        Data.MSFPAcq.EXC = num_EXC;
        Data.MSFPAcq.ROI = num_ROI;
        
        
        % Previously, because the Doric software names the first ROI as
        % ROI0. Mine is 0,2,3,4,5,6,7,8,9,10,10,11. But after batches after
        % 20221006, the problem seems to be fixed
        if num_CAM == 1
            for i_column = 1:num_var
                aux = strjoin({names1{i_column}(1:min(length(names1{i_column}),6)),'_',names2{i_column}});
                var_names{i_column} = aux(~isspace(aux));             
            end
        elseif num_CAM == 2
            for i_column = 1:num_var
            
                if i_column == 8 | i_column == 9 | i_column == 10
                    aux = strjoin({names1{i_column}(1:min(length(names1{i_column}),6)),'_',names2{i_column}});
                    var_names{i_column} = strrep(aux(~isspace(aux)),'0','1');
                   % if you have 10 fibers
                elseif i_column == 38 | i_column == 39 | i_column == 40
                    aux = strjoin({names1{i_column}(1:min(length(names1{i_column}),6)),'_',names2{i_column}});
                    var_names{i_column} = strrep(aux(~isspace(aux)),'10','11');
                elseif i_column == 41 | i_column == 42 | i_column == 43
                    aux = strjoin({names1{i_column}(1:min(length(names1{i_column}),6)),'_',names2{i_column}});
                    var_names{i_column} = strrep(aux(~isspace(aux)),'11','12');
                    
                else
                    aux = strjoin({names1{i_column}(1:min(length(names1{i_column}),6)),'_',names2{i_column}});
                    var_names{i_column} = aux(~isspace(aux));
                
                end
            end

        end
        
        if num_CAM == 2

            if DoricStudioVersion == '5.4.1.23'  % Header labels are swapped in this version of the Doric Studio software, so we are correcting this here:
                fprintf('Variable names before correction:\n')
                disp(var_names(1:8));
        
                indCAM2 = strfind(var_names, 'CAM#2');
                IndexCAM2 = find(not(cellfun('isempty',indCAM2)));
                aux = var_names{IndexCAM2};
                indEXC1 = strfind(var_names, 'EXC#1');
                IndexEXC1 = find(not(cellfun('isempty',indEXC1)));
                indEXC2 = strfind(var_names, 'EXC#2');
                IndexEXC2 = find(not(cellfun('isempty',indEXC2)));
                indEXC3 = strfind(var_names, 'EXC#3');
                IndexEXC3 = find(not(cellfun('isempty',indEXC3)));
        
                var_names{IndexCAM2} = var_names{IndexEXC1};  % CAM2 = EXC1
                var_names{IndexEXC1} = var_names{IndexEXC2};  % EXC1 = EXC2
                var_names{IndexEXC2} = var_names{IndexEXC3};  % EXC2 = EXC3
                var_names{IndexEXC3} = aux;                   % EXC3 = CAM2
                
                fprintf('Variable names AFTER correction:\n')
                disp(var_names(1:8));
            end
        
        end

        % Read the rest of the table
        doric_table = readtable([doric_data_folder doric_filename{1}],'HeaderLines',2);  % skips the first two rows of data
        time_acq = doric_table{:,1};
    
        % Create two cells: one for time, one for the variable values. Each column
        % is one variable
        var_time = cell(1,num_var);
        var_values = cell(1,num_var);
        if plotauxfig
            figure; hold on
        end
        for i_column = 2:num_var
            aux = doric_table{:,i_column};
            var_time{i_column} = time_acq(~isnan(aux));   % not all variables are generated/measured at the same time points (time_acq includes all possible time points in which something is generated/measured)
            var_values{i_column} = aux(~isnan(aux));
    
            if plotauxfig && i_column <= 10
                subplot(12,1,i_column-1); hold on
        %         plot(var_time{i_column}(9000:19000),var_values{i_column}(9000:19000))
                plot(var_time{i_column},var_values{i_column})
                ylabel(var_names{i_column})
                if i_column == length(names1)
                    xlabel('Time (s)')
                    ylabel(var_names{i_column})
                end  
            end
        end
    
        clear('doric_table')
    

    
        %% Plot the cameras' acquisition and the corresponding excitations
    
        if plotauxfig && num_EXC == 2
            indC1 = strfind(var_names, 'CAM#1');
            IndexC1 = find(not(cellfun('isempty',indC1)));
            
            indE1 = strfind(var_names, 'EXC#1');
            IndexE1 = find(not(cellfun('isempty',indE1)));
            indE2 = strfind(var_names, 'EXC#2');
            IndexE2 = find(not(cellfun('isempty',indE2)));
            
    
            figure;
            subplot(1,1,1); hold on % Camera 1 and its excitations (1 and 2):
            plot(var_time{IndexC1},var_values{IndexC1},'o')
            plot(var_time{IndexE1},var_values{IndexE1},'color',[221 85 25]/255)
            plot(var_time{IndexE2},var_values{IndexE2},'--g')
            ylim([-0.2 1.2])
            xlim([1096 1097])
            legend('CAM1','415nm','470nm')

        elseif plotauxfig && num_EXC == 3
            indC1 = strfind(var_names, 'CAM#1');
            IndexC1 = find(not(cellfun('isempty',indC1)));
            indC2 = strfind(var_names, 'CAM#2');
            IndexC2 = find(not(cellfun('isempty',indC2)));
            indE1 = strfind(var_names, 'EXC#1');
            IndexE1 = find(not(cellfun('isempty',indE1)));
            indE2 = strfind(var_names, 'EXC#2');
            IndexE2 = find(not(cellfun('isempty',indE2)));
            indE3 = strfind(var_names, 'EXC#3');
            IndexE3 = find(not(cellfun('isempty',indE3)));
    
            figure;
            subplot(2,1,1); hold on % Camera 1 and its excitations (1 and 2):
            plot(var_time{IndexC1},var_values{IndexC1},'o')
            plot(var_time{IndexE1},var_values{IndexE1},'color',[221 85 25]/255)
            plot(var_time{IndexE2},var_values{IndexE2},'--g')
            ylim([-0.2 1.2])
            xlim([1096 1097])
            legend('CAM1','415nm','470nm')
            subplot(2,1,2); hold on % Camera 2 and its excitation (3):
            plot(var_time{IndexC2},var_values{IndexC2},'o')
            plot(var_time{IndexE3},var_values{IndexE3},'r')
            ylim([-0.2 1.2])
            xlim([1096 1097])
            legend('CAM2','568nm')
    
            figure; hold on
            plot(var_time{IndexC1},var_values{IndexC1},'o')
            plot(var_time{IndexE1},var_values{IndexE1},'color',[221 85 25]/255)
            plot(var_time{IndexE2},var_values{IndexE2},'g')
            plot(var_time{IndexC2},var_values{IndexC2},'o')
            plot(var_time{IndexE3},var_values{IndexE3},'m')
            ylim([-0.2 1.2])
            xlim([1096 1097])
            legend('CAM1','415nm','470nm','CAM2','568nm')
    
    
        end
    
        %% Check time between data points in each fluotrace
    
        % Check reliability of acquisition and save frame rates
        Data.MSFPAcq.FrameRate = [];
        if plotauxfig; figure; hold on; end
        for i_exc = 1:num_EXC
            ind = strfind(var_names, ['EXC' num2str(i_exc)]);
            Index = find(not(cellfun('isempty',ind)));
            Index = Index(1);   % look only at the times of acquisition of the 1st ROI (all other ROIs are acquired at the same time)
    
            checkFR = var_time{Index}(2:end)-var_time{Index}(1:end-1);
    
            if plotauxfig
                subplot(1,num_EXC,i_exc); histogram(checkFR)
                title ('Time between frames')
                legend(['Frame Rate:' num2str(1/mode(checkFR))])
                title (['EXC' num2str(i_exc)])
            end
            Data.MSFPAcq.FrameRate = [Data.MSFPAcq.FrameRate; 1/mode(checkFR)];
        end
    
        %% Plot all acquisitions of each ROI
    
        figure
        for i_ROI = 1:num_ROI
            ind = strfind(var_names, ['ROI' num2str(i_ROI-1)]);  
            Index = find(not(cellfun('isempty',ind)));
    
            subplot(ceil(num_ROI/6),6,i_ROI); hold on
            for i_EXC = 1:num_EXC
                plot(var_time{Index(i_EXC)},var_values{Index(i_EXC)}-mean(var_values{Index(i_EXC)}))
            end
            ylabel(['ROI' num2str(i_ROI-1)])
            if i_ROI >= ceil(num_ROI/6)*5
                xlabel('Time (s)')
            end
            
            if i_ROI == 1
                if num_CAM == 1
                    legend('Exc1','Exc2')
                elseif num_CAM == 2
                    legend('Exc1','Exc2','Exc3')
                end
            end
        end
    
        %% Find the synch pulses (to align to behavior)
    
        ind = strfind(var_names, 'DI/O');
        Index = find(not(cellfun('isempty',ind)));
        if numel(Index) > 1; Index = Index(1); end
        DIO = var_values{Index};
        if strcmp(filename{1,1},'FgDA_01_Selina_C5D5R3E5R3_20220214_cond6.mat')
            DIO(800000:end) = 0;
        end
        DIO_time = var_time{Index};
    
        x = DIO/max(DIO);
        % detect trial start
        topThres = 0.8;
        x(1) = topThres+0.1;
        a = find(x > topThres);
        b = [0; diff(a)]; 
        TTLsynch1 = a(b > 1);
        TTLsynchTimesUp = DIO_time(TTLsynch1);
        a = find(x < 0.1);
        b = [0; diff(a)]; 
        TTLsynch2 = a(b > 1);
        TTLsynchTimesDown = DIO_time(TTLsynch2);
        
        num_TTL = numel(TTLsynch1);
    
        figure; hold on
        plot(DIO_time,DIO)
        plot(TTLsynchTimesUp,0.8,'or')
        plot(TTLsynchTimesDown,0.8,'+y')
    
        %% Make a cell that separates the ROI values by trial
    
        ind = strfind(var_names, 'ROI');
        Index = find(not(cellfun('isempty',ind)));
        Index = Index(1):Index(1)+num_EXC-1;   % look only at the times of acquisition of the 1st ROI (all other ROIs are acquired at the same time)
    
        time_acq_ROIsFluos = cell(1,num_EXC);
        time_trials_ROIsFluos = cell(1,num_TTL);
        fluo_trials = cell(1,num_TTL);
        
    %     % Each excitation will have its own acquisition time resulting in 3 variables for each ROI
    %     % To simplify things, I will create a fluo time vector for each trial
    %     % and interpolate the fluo values of each EXC
    % 
    %     for i_trial = 1:num_TTL
            
        
        for i_exc = 1:num_EXC                     
            
            % Build a camera time vector that is scaled to compensate for time
            % drift between the CAM clock and the DIO(BFPD) clock:
            clocks_scalingfactor = DIO_time(end)/var_time{Index(i_exc)}(end);
            
            time_acq_ROIsFluos{i_exc} = clocks_scalingfactor*var_time{Index(i_exc)};
    
            for i_trial = 1: min(length(TTLsynch2),159)
                
                % Instead of taking the time vector for each fluo, I am going
                % to interpolate on a defined time vector (based on the DIO) for all fluos
                
                time_interp = TTLsynchTimesUp(i_trial):1/Data.MSFPAcq.FrameRate(i_exc):TTLsynchTimesDown(i_trial);
    
                p1 = find(time_acq_ROIsFluos{i_exc}-TTLsynchTimesUp(i_trial)>0,1,'first');
                
                if numel(TTLsynchTimesDown)>= i_trial
                    p2 = find(time_acq_ROIsFluos{i_exc}-TTLsynchTimesDown(i_trial)<0,1,'last');
                else
                    p2 = find(time_acq_ROIsFluos{i_exc}-DIO_time(end)<0,1,'last');
                end

                
                for i_roi = 1:num_ROI
                    index_roiexc = strfind(var_names, ['ROI' num2str(i_roi-1)]); 
                    Index_roi = find(not(cellfun('isempty',index_roiexc)));
                    if time_acq_ROIsFluos{i_exc}(p1) > time_acq_ROIsFluos{i_exc}(p1+1)
                        time_acq_ROIsFluos{i_exc}(p1) = time_acq_ROIsFluos{i_exc}(p1+1)-0.04;
                    end
                    fluo_interp = interp1(time_acq_ROIsFluos{i_exc}(p1:p2),var_values{Index_roi(i_exc)}(p1:p2),time_interp);
                    time_trials_ROIsFluos{i_trial}{i_exc} = time_interp;
                    fluo_trials{i_trial}{i_exc}(:,i_roi) = fluo_interp;    % this vector will be used to match with behavior data
                end
            end
        end
    else
        num_ROI = 10;
        fluodataAvailable = 0;
        fakeEfctFR = 20;  % fake effective frame rate
        TTLsynch1 = [];
    end
    
    %% Load corresponding behavioral session and put it together with MSFP
    
    % beh_data_folder = [''];
    % oldfolder = cd(beh_data_folder);
    % session = dir(['*' day '*']);
    % load([session.folder '/' session.name]);   % this loads a strcuture called 'SessionData'
    % session.name

    beh_data_folder = [strcat(filedir,'/bpod/')];
    beh_filename = [filename{1}];
    load([beh_data_folder beh_filename]);
    
    
    Data.Session.TrialTypes = SessionData.TrialTypes;
    
    Data.Session.Rewards = SessionData.RewardDelivered; 

    if length(TTLsynch1) ~= SessionData.nTrials &&...
            length(TTLsynch1) ~= SessionData.nTrials+1 &&...
            length(TTLsynch1) ~= SessionData.nTrials-1 && fluodataAvailable
        disp('Unmatched number of trials')
    else
        num_trials = min(SessionData.nTrials,length(TTLsynch1));
        for i_trial = 1:num_trials
            % Behavioral data
    %         Data.RotEncoder(i_trial) = SessionData.RotaryEncoderData(i_trial);
            if ~isfield(SessionData.RawEvents.Trial{i_trial}.Events,'Port1Out')
                SessionData.RawEvents.Trial{i_trial}.Events.Port1Out = [];
            end
            if ~isfield(SessionData.RawEvents.Trial{i_trial}.Events,'Port1In')
                SessionData.RawEvents.Trial{i_trial}.Events.Port1In = [];
            end
            Data.Events(i_trial) = SessionData.RawEvents.Trial{i_trial}.Events;
            Data.States(i_trial) = SessionData.RawEvents.Trial{i_trial}.States;
            if fluodataAvailable
                % Fiber-photometry data
                Data.MSFP(i_trial).Time = time_trials_ROIsFluos{i_trial};
                Data.MSFP(i_trial).Fluo = fluo_trials{i_trial};
            else
                % In case I forgot to save the MSFP data (pressed 'Live'
                % instead of 'Record' in the Doric Software
                Data.MSFP(i_trial).Time = 0:1/fakeEfctFR:SessionData.RawEvents.Trial{i_trial}.States.ITI(end)-1/fakeEfctFR;
                Data.MSFP(i_trial).Fluo = nan(size(Data.MSFP(i_trial).Time));
            end
        end
    end
    
    %% Saving
    
    % Data has fields: FrameRate, RotEncoderTime, RotEncoderPositions
    % Events, States, MSFP
    
    foldername_tosave = [filedir '/processed/'];
    filename_tosave = [filename{1}];
    save([foldername_tosave filename_tosave],'-struct','Data');  % make sure to use '-struct' as argument to save each field of 'Data' independently

end
