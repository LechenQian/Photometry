function get_autocorrelation_rawMSFP(mouse,day,DoricStudioVersion)

dbstop if error

% e.g. mouse: 'SM79'
% e.g. day: '20190404'
plotauxfig = 1;


%% Load the data saved by the Doric Studio software (from the text file .csv)
% I will load the .csv file and look at the output of the ROI in the
% different cameras/acquisition and test the movement artifact correction
% method that is more appropriate

doric_data_folder = '/home/sara/Documents/DATA/MultiSiteFiberPhotometry/Data_Doric/';
doric_filename = [mouse '_' day '_0000.csv'];

% Get the names of the variables that describe each column (two first lines
% of the file)
% Note: the labels might not be always the exact same in all .csv files,
% depending on the DI/O channel was selected to be included, for instance

cd(doric_data_folder);     % Returns all the files and folders in the directory                  
doric_file = dir(['*' mouse '_' day '*']);

if ~isempty(doric_file)

    fluodataAvailable = 1;
    fprintf('Processing session %s\n',doric_file.name);
    
    % Load headers on the first row
    fid = fopen([doric_data_folder doric_filename], 'r');
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
    for i_column = 1:num_var
        aux = strjoin({names1{i_column}(1:min(length(names1{i_column}),6)),'_',names2{i_column}});
        var_names{i_column} = aux(~isspace(aux));
    end

    if strcmp(DoricStudioVersion, '5.4.1.23')  % Header labels are swapped in this version of the Doric Studio software, so we are correcting this here:
        fprintf('Variable names before correction:\n')
        disp(var_names(1:8));
        if strcmp(mouse,'SM144') && strcmp(day,'20211210')  % on this day I saved 2 DIO inputs by mistake and I this got the saving of headers wrong. Here's the fix:
            
            indDIO1 = strfind(var_names, 'DI/O#_');
            IndexDIO1 = find(not(cellfun('isempty',indDIO1)));            
            aux = var_names{IndexDIO1(1)};
            indCAM1 = strfind(var_names, 'CAM#1');
            IndexCAM1 = find(not(cellfun('isempty',indCAM1)));
            indCAM2 = strfind(var_names, 'CAM#2');
            IndexCAM2 = find(not(cellfun('isempty',indCAM2)));
            indEXC1 = strfind(var_names, 'EXC#1');
            IndexEXC1 = find(not(cellfun('isempty',indEXC1)));
            indEXC2 = strfind(var_names, 'EXC#2');
            IndexEXC2 = find(not(cellfun('isempty',indEXC2)));
            indEXC3 = strfind(var_names, 'EXC#3');
            IndexEXC3 = find(not(cellfun('isempty',indEXC3)));
            
            aux2 = var_names{IndexCAM2};
            var_names{IndexDIO1(1)} = var_names{IndexCAM1};   % DIO#1 = CAM1
            var_names{IndexCAM1} = var_names{IndexEXC1};      % CAM1 = EXC1
            var_names{IndexCAM2} = var_names{IndexEXC2};      % CAM2 = EXC2
            var_names{IndexEXC1} = aux2;                      % EXC1 = CAM2
            var_names{IndexEXC2} = var_names{IndexEXC3};      % EXC2 = EXC3
            var_names{IndexEXC3} = aux;                       % EXC3 = DIO#1     
        else
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
        end
        fprintf('Variable names AFTER correction:\n')
        disp(var_names(1:8));
    end

    % Read the rest of the table
    doric_table = readtable([doric_data_folder doric_filename],'HeaderLines',2);  % skips the first two rows of data
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

    %% Plot the cameras' acquisition and the corresponding excitations

    if plotauxfig && num_EXC == 3
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
        legend('CAM2','40nm')

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
        ind = strfind(var_names, ['ROI' num2str(i_ROI-1)]);  % Because the Doric software names the first ROI as ROI0
        Index = find(not(cellfun('isempty',ind)));

        subplot(ceil(num_ROI/6),6,i_ROI); hold on
        for i_EXC = 1:num_EXC
            plot(var_time{Index(i_EXC)},var_values{Index(i_EXC)}-mean(var_values{Index(i_EXC)}))
        end
        ylabel(['ROI' num2str(i_ROI)])
        if i_ROI >= ceil(num_ROI/6)*5
            xlabel('Time (s)')
        end
        if i_ROI == 1
            legend('Exc1','Exc2','Exc3')
        end
    end
    
    % Figure for lab meeting
    if plotauxfig
        MCoorTBL = table_MouseImplantFibers_MSFP;
        fiberIDs = table2cell(MCoorTBL(:,'BrainRegion'));
        fiberIDs = cellfun(@(x) char(x),fiberIDs,'UniformOutput',false);
        figure('Position',[675 373 1066 601])
        set(gcf , 'Color', 'w'); hold on
        for i_ROI = 1:num_ROI
            ind = strfind(var_names, ['ROI' num2str(i_ROI-1)]);  % Because the Doric software names the first ROI as ROI0
            Index = find(not(cellfun('isempty',ind)));

            subplot(1,3,1); hold on      
            plot(var_time{Index(1)},var_values{Index(1)}-mean(var_values{Index(1)})+i_ROI*1200,'linewidth',1)
            text(1249, var_values{Index(2)}(1)-mean(var_values{Index(2)})+i_ROI*1200,fiberIDs(i_ROI))
            if i_ROI == num_ROI
                xlim([1255 1280])
                set(gca,'Tickdir','out', 'Xcolor', 'w', 'Ycolor', 'w')
                title('Isosbestic')
            end
            subplot(1,3,2); hold on      
            plot(var_time{Index(2)},var_values{Index(2)}-mean(var_values{Index(2)})+i_ROI*1200,'linewidth',1)
            if i_ROI == num_ROI
                xlim([1255 1280])
                set(gca,'Tickdir','out', 'Xcolor', 'w', 'Ycolor', 'w')
                title('GCaMP')
            end
            subplot(1,3,3); hold on      
            plot(var_time{Index(3)},var_values{Index(3)}-mean(var_values{Index(3)})+i_ROI*1200,'linewidth',1)
            if i_ROI == num_ROI
                xlim([1255 1280])
                set(gca,'Tickdir','out', 'Xcolor', 'w', 'Ycolor', 'w')
                plot([1255 1255],[18080 18080+1000],'color','k','linewidth',1)
                text(1249, 18080+500,'1000au')
                plot([1256 1261],[18080+1000 18080+1000],'color','k','linewidth',1)
                text(1257, 18080+1500,'5s')
                title('tdTomato')
            end
        end       
        
        % Figure to look at autocorrelation of signals across different
        % areas:
        % Calculate deltaF/F based on a moving mean
        
        
        ftraces = figure('Position',[475 1 1562 970]);
        set(gcf , 'Color', 'w'); hold on
        fxcorr = figure('Position',[675 373 1066 601]);
        set(gcf , 'Color', 'w'); hold on
        FR = round(max(Data.MSFPAcq.FrameRate));
        autoCorr = nan(num_EXC,num_ROI,FR*5*2+1);
        autoCorr_t = nan(num_EXC,num_ROI,FR*5*2+1);
        for i_ROI = 1:num_ROI
            ind = strfind(var_names, ['ROI' num2str(i_ROI-1)]);  % Because the Doric software names the first ROI as ROI0
            Index = find(not(cellfun('isempty',ind)));
            
            prctl = 6;
            
            for i_exc = 2 %1:num_EXC
                fluo = var_values{Index(i_exc)};
                fluot_time = var_time{Index(i_exc)};
                F0_fluo = running_percentile(fluo, FR*40, prctl);  % Calculate F0 (40s running 6th percentile), of both neuropil corrected and raw fluo
                DFoF = (fluo - F0_fluo)./F0_fluo;
                DFoF_norm = DFoF/max(DFoF);
                
                figure(ftraces)
%                 subplot(1,3,i_exc); hold on
                plot(fluot_time,DFoF_norm+i_ROI*1.1,'linewidth',1)
%                 xlim([1255 1280])
                xlim([1180 1280])
                if i_exc == 2
                    text(1175, DFoF_norm(1)+i_ROI*1.1,fiberIDs(i_ROI))
                end
                if i_ROI == num_ROI
                    set(gca,'Tickdir','out', 'Xcolor', 'w', 'Ycolor', 'w')
                    switch i_exc
                        case 1
                            title('Isosbestic')
                        case 3
                            title('tdTomato')
                        case 2
                            title('GCaMP')
                            set(gca,'Tickdir','out', 'Xcolor', 'w', 'Ycolor', 'w')
                            plot([1180 1180],[17 17+1],'color','k','linewidth',1)
                            text(1181, 17+0.5,'1 norm\DeltaF/F')
                            plot([1180 1180+10],[18 18],'color','k','linewidth',1)
                            text(1183, 18+0.5,'10s')
                    end
                end
                
                %% Get autocorrelogram
                switch i_exc
                    case 1
                        cor = [0.7 0.7 0.7];
                    case 2
                        cor = [0.1 0.8 0.1];
                    case 3
                        cor = [1 0.7 0.7];
                end
                [autocorr, lags] = xcorr(DFoF_norm,FR*5,'coeff');
                autoCorr(i_exc,i_ROI,:) = autocorr;
                autoCorr_t(i_exc,i_ROI,:) = lags/FR;
                figure(fxcorr)
                subplot(3,5,i_ROI); hold on
                plot(lags/FR,autocorr,'color',cor,'linewidth',2)
                xlabel('Lag (s)')
                ylabel('Autocorrelation')
                ylim([0.4 1])
                title(fiberIDs(i_ROI))
            end 
        end
    end
    
    %% Plot autocorrelation traces of all areas together
    FiberColors = brewermap(num_ROI,'*Spectral');
    exp_fun = fittype('a+(1-a)*exp(-(1/c)*x)');
    fit_start = FR*5+1;
    fit_end = FR*5+1+FR*2;
    
    for i_exc = 2 %1:num_EXC
        figure('Name',['Excitation' num2str(i_exc)],'Position',[675 521 1240 450]);
        subplot(1,3,1); hold on   
        h1 = plot(squeeze(autoCorr_t(i_exc,:,:))',squeeze(autoCorr(i_exc,:,:))','linewidth',3);
        xlim([0 1.6])
        set(h1,{'color'},num2cell(FiberColors,2))
        set(gca,'TickDir','out')
        xlabel('Lag(s)')
        ylabel('Autocorrelation')
        legend(fiberIDs,'NumColumns',3)

        % Fit and exponential to the decay of the autocorrelation:
        subplot(1,3,3); hold on
        plot([4.5 4.5],[0 0.6],'k')
        plot([7.5 7.5],[0 0.6],'k')
        plot([10.5 10.5],[0 0.6],'k')
        for i_ROI = 1:num_ROI
            xfit = squeeze(autoCorr_t(i_exc,i_ROI,fit_start:fit_end));
            yfit = squeeze(autoCorr(i_exc,i_ROI,fit_start:fit_end));
            yfit_norm = (yfit-min(yfit))/(max(yfit)-min(yfit));
            try
                [fit_exp,~,~] = fit(xfit,yfit_norm,exp_fun,'StartPoint',[0.2 100]);
                slopes_exp(i_exc,i_ROI) = fit_exp.c;
                fitted_exp(i_exc,i_ROI,:) = fit_exp(squeeze(autoCorr_t(i_exc,i_ROI,fit_start:fit_end)));
            catch
                slopes_exp(i_exc,i_ROI) = nan;
                fitted_exp(i_exc,i_ROI,:) = nan(size(squeeze(autoCorr_t(i_exc,i_ROI,fit_start:fit_end))));
            end
            subplot(1,3,2); hold on
            plot(xfit,yfit_norm,'color',FiberColors(i_ROI,:),'linewidth',2);
            plot(squeeze(autoCorr_t(i_exc,i_ROI,fit_start:fit_end)),squeeze(fitted_exp(i_exc,i_ROI,:)),':k','linewidth',0.5,'HandleVisibility','off');
            set(gca,'TickDir','out')
            ylabel('Autocorrelation (normalized)')
            xlabel('Lag(s)')
        
            subplot(1,3,3); hold on
            bar(i_ROI,slopes_exp(i_exc,i_ROI),'FaceColor',FiberColors(i_ROI,:),'EdgeColor','k');
        end
        xticks(1:num_ROI);
        xtickangle(45)
        ylabel('Time constant(s)')
        set(gca,'TickDir','out','xticklabels',fiberIDs)
        
    end      
end
end