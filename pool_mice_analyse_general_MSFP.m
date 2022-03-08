function pool_mice_analyse_general_MSFP(protocol,taskPhase,alignTo,smooth_time,smooth_acrosstrials,sinal_allplots,sinal_filtLicking,mov_correction)
dbstop if error
rng(2020)
% Explanation of the inputs required to run this function:

% protocol: 'DRL6Rv' or 'Revers' or 'DRL6Od' or 'VarMag' or 'DRL3Od', 'WDMSFP'  (see 'Get protocol specific info' section below)
% taskPhase: '0' initial learning; '1' expert (before any putative
% reversal), '2' reversal learning day; '3' expert after reversal learning
% alignTo: 'CS' (or 'US')
% smooth_time: to smooth in time
% smooth_acrosstrials: yes or no
% sinal_allplots: 1 or 0 to do all plots or not (including auxiliary to see all cells)
% sinal_filtLicking: 1 for not accounting for trials after the mouse has stopped licking
% mov_correction: whether to perform movement correction or not

%% Target areas of each fiber:

MCoorTBL = table_MouseImplantFibers_MSFP;
fiberIDs = table2cell(MCoorTBL(:,'BrainRegion'));
fiberIDs = cellfun(@(x) char(x),fiberIDs,'UniformOutput',false);

%% Get protocol specific info and data (I should separate this part and put it in another general function to call the entire analysis from)

[prt_colors, prt_conditions, num_CSTT, ~, ~] = get_protocol_info(protocol);
fluotype = 'FiberPhotometry';

switch protocol
    case 'WDMSFP'
        sinal_filtLicking = 0;
        taskPhase = 1;
        mouse = {'SM143','SM143','SM144','SM145'};
        day = {'20211117','20211130','20211130','20211130'};
    case 'Es163O'
        sinal_filtLicking = 0;
        taskPhase = 1;
        switch taskPhase
            case 0
                mouse = {'SM143','SM143','SM145','SM145','SM144'};
                day = {'20211201','20211202','20211201','20211202','20211202'};
            case 1
        end
    case 'Eshl16'
        sinal_filtLicking = 1;
        taskPhase = 1;
        session = {'SM143_20211210','SM144_20211210','SM145_20211210'};
end
% Instead of doing this, I might want to change the code below to always
% work in terms of session ? Decide what to do later
if ~exist('mouse') && exist('session')
    for i_session = 1:numel(session)
        aux = strsplit(session{i_session},'_');
        mouse{i_session}  = aux{1};
        day{i_session} = aux{2};
    end
end
    
%% Load the data (all trials in a session separated by trial type)
num_mice = length(mouse);
num_TT = nan(1,num_mice);
num_cells = nan(1,num_mice);
licks_all = cell(1,num_mice);
fluo_all = cell(1,num_mice);
revtrial = cell(1,num_mice);
TT_allRsizes = cell(1,num_mice);
Rsize_pos = cell(1,num_mice);
CSin = cell(1,num_mice);
CSout = cell(1,num_mice);
TraceEnd = cell(1,num_mice);
SR_fluo = cell(1,num_mice);

SR_lick = 1000;
t1 = tic;
parfor i_mouse = 1:num_mice
    [licks_all{i_mouse}, time_licks{i_mouse}, fluo_all{i_mouse}, ...
        time_fluo{i_mouse}, revtrial{i_mouse}, TT_allRsizes{i_mouse},...
        Rsize_pos{i_mouse}, CSin{i_mouse}, CSout{i_mouse}, ...
        TraceEnd{i_mouse}, SR_fluo{i_mouse},num_EXC{i_mouse}, ...
        RotEnc_all{i_mouse}] = load_session_data_MSFP(protocol,...
        mouse{i_mouse}, day{i_mouse}, alignTo, sinal_filtLicking);
    % dim of fluo_all: time,trials,rois,wavelength excitations
    % revtrial contains the first trial with reversed outcome for each TT with a CS

    num_TT(i_mouse) = length(licks_all{i_mouse});          % total number of trial types for each mouse
    num_cells(i_mouse) = size(fluo_all{i_mouse}{1},3);     % total number of cells for each mouse
end
if taskPhase ~= 2   % To use the same code for reversal and other protocols (in this case set reversal trial to 1)
    revtrial = cellfun(@(x) ones(size(x)),revtrial, 'UniformOutput',false);
end
toc(t1)

time_licks = time_licks{1};
time_fluo = time_fluo{1};
time_rotenc = time_fluo;

%% Perform movement artifact correction

if num_EXC{1} == 3
    funcChannel = 2;
    useITIonly = 1;
    if mov_correction
        if useITIonly 
            % Get the indices for the regions of ITI and use that for
            % movement correction
            CSstartT = 0; CSendT = 1;
            [~, ind_CS, ~, ~, ind_US_exp_early] = get_indicesEvents(time_fluo,fluotype,CSin,CSstartT,CSendT,TraceEnd);
            ind_touse = cat(2,1:ind_CS-10,ind_US_exp_early(end)+1:numel(time_fluo));
%             ind_touse = ind_US_exp_early(1)+1:numel(time_fluo);
        else 
            % Use the entire trial duration for movement correction
            ind_touse = 1:numel(time_fluo);
        end       
        fluo_all_reduc = get_movcorrected_fluo(fluo_all,ind_touse,funcChannel,mouse,day,fiberIDs);
    else
        % Change the shape of the cells, to include only the functional channel
        % (so that all the functions written for 2P also work here
        
        % fluo_all{i_mouse} has one cell per trial type, each with dimensions
        % time, trials, ROIs, wavelengths. I will make of the dimensions time,
        % trials, ROIs:
        fluo_all_reduc = cell(size(fluo_all));
        for i_mouse = 1:num_mice
            fluo_all_reduc{i_mouse} = cellfun(@(x) squeeze(x(:,:,:,funcChannel)),fluo_all{i_mouse}, 'UniformOutput',false);
        end
    end
    fluo_all = fluo_all_reduc;
    clear fluo_all_reduc
end

%% Smooth fluo data (in time)

fluo_all_filt = cell(size(fluo_all));
if smooth_time
    % filter along time
    N = SR_fluo{1};   % 1s
%     alpha = 5;        % (30-1)/(2*5)
    alpha = 15;
    g1 = gausswin(N,alpha)/sum(gausswin(N,alpha));
    
    % apply filtering to all TT of each mouse simultaneously
    for i_mouse = 1:num_mice
        fluo_all_filt{i_mouse} = cellfun(@(x) convn(x,g1,'same'),fluo_all{i_mouse}, 'UniformOutput',false);
    end
    fluo_all = fluo_all_filt;
    clear fluo_all_filt
end


%% Calculate the mean response of each ROI to each TT. Plot this mean and trials of each cell in the session

flicks = figure('name','Mean licks traces','Position',[100 1000 1000 200]); hold on
frotenc = figure('name','Rotary encoder traces','Position',[100 1000 1000 650]); hold on
rotencvar = {'Distance (m)','Velocity (m/s)','Acceleration (m/s^2)'};
for i_mouse = 1:num_mice

    lick_mean_timecourse{i_mouse} = cellfun(@(x) squeeze(nanmean(x,2)),licks_all{i_mouse}, 'UniformOutput',false);  % dim: time,cell
    lick_std_timecourse{i_mouse} = cellfun(@(x) squeeze(nanstd(x,0,2)),licks_all{i_mouse}, 'UniformOutput',false);
    
    fluo_mean_timecourse{i_mouse} = cellfun(@(x) squeeze(nanmean(x,2)),fluo_all{i_mouse}, 'UniformOutput',false);   % dim: time,ROI,fluo wavelength
    fluo_std_timecourse{i_mouse} = cellfun(@(x) squeeze(nanstd(x,0,2)),fluo_all{i_mouse}, 'UniformOutput',false);
    
    rotenc_mean_timecourse{i_mouse} = cellfun(@(x) squeeze(nanmean(x,2)),RotEnc_all{i_mouse}, 'UniformOutput',false);   % dim: time,rot enc info: position/velicity/acceleration
    rotenc_std_timecourse{i_mouse} = cellfun(@(x) squeeze(nanstd(x,0,2)),RotEnc_all{i_mouse}, 'UniformOutput',false);
        
    % plot mean licking pattern of each mouse (all TTs together)
    smoothlicks = gausswin(1000,3)/sum(gausswin(1000,3));
    figure(flicks);
    subplot(1,num_mice,i_mouse) ; hold on
    title([mouse{i_mouse} '_' day{i_mouse}])
    line([CSin{i_mouse} CSin{i_mouse}],[0 15],'color','k')
    line([TraceEnd{i_mouse} TraceEnd{i_mouse}],[0 15],'color','k')
    for i_TT = 1:num_CSTT
        meantoplot = conv(lick_mean_timecourse{i_mouse}{i_TT},smoothlicks,'same');
        semtoplot = conv(lick_std_timecourse{i_mouse}{i_TT},smoothlicks,'same')/sqrt(size(licks_all{i_mouse}{i_TT},2));
        shadedErrorBar(time_licks,meantoplot,semtoplot,{'color',prt_colors(i_TT,:),'linewidth',2},0.2)
        xlim([-1 6])
        if i_mouse == 1
            ylabel('Lick rate (s^{-1})')
            xlabel('Time-Cue (s)')
        end
    end
    
    % plot mean Rotary encoder pattern of each mouse (all TTs together)
    for i_rotencvar = 1:size(rotencvar,2)
        figure(frotenc);
        subplot(3,num_mice,i_mouse+(i_rotencvar-1)*num_mice) ; hold on
        title([mouse{i_mouse} '_' day{i_mouse}])
        line([CSin{i_mouse} CSin{i_mouse}],[0 15],'color','k')
        line([TraceEnd{i_mouse} TraceEnd{i_mouse}],[0 15],'color','k')
        for i_TT = 1:num_CSTT
    %         meantoplot = conv(lick_mean_timecourse{i_mouse}{i_TT},smoothlicks,'same');
    %         semtoplot = conv(lick_std_timecourse{i_mouse}{i_TT},smoothlicks,'same')/sqrt(size(licks_all{i_mouse}{i_TT},2));
            meantoplot = rotenc_mean_timecourse{i_mouse}{i_TT}(:,i_rotencvar);
            semtoplot = rotenc_std_timecourse{i_mouse}{i_TT}(:,i_rotencvar);
            shadedErrorBar(time_rotenc,meantoplot,semtoplot,{'color',prt_colors(i_TT,:),'linewidth',2},0.2)
            xlim([-1 6])
            if i_mouse == 1
                ylabel(rotencvar(i_rotencvar))
                xlabel('Time-Cue (s)')
            end
        end
    end

    if sinal_allplots
        % All trials in a session (actually I set a limint for initial and reversal learning)
        plot_pcolorplot_percell(taskPhase,mouse{i_mouse},day{i_mouse},licks_all{i_mouse},time_licks,lick_mean_timecourse{i_mouse},lick_std_timecourse{i_mouse},prt_colors,prt_conditions,CSin{i_mouse},CSout{i_mouse},TraceEnd{i_mouse},revtrial{i_mouse},SR_lick,0,num_CSTT,fiberIDs);
        plot_pcolorplot_percell(taskPhase,mouse{i_mouse},day{i_mouse},fluo_all{i_mouse},time_fluo,fluo_mean_timecourse{i_mouse},fluo_std_timecourse{i_mouse},prt_colors,prt_conditions,CSin{i_mouse},CSout{i_mouse},TraceEnd{i_mouse},revtrial{i_mouse},SR_fluo{i_mouse},1,num_CSTT,fiberIDs);
        plot_pcolorplot_percell(taskPhase,mouse{i_mouse},day{i_mouse},RotEnc_all{i_mouse},time_rotenc,rotenc_mean_timecourse{i_mouse},rotenc_std_timecourse{i_mouse},prt_colors,prt_conditions,CSin{i_mouse},CSout{i_mouse},TraceEnd{i_mouse},revtrial{i_mouse},SR_fluo{i_mouse},1,num_CSTT,fiberIDs);
        
        % Only trials where learning is rapidly happening
        %         if strcmp(protocol,'Revers') || strcmp(protocol,'DRL6Rv')
        %             plot_revers_closeup_percell(fluo_all{i_mouse},time_fluo,revtrial{i_mouse},CSin{i_mouse},TraceEnd{i_mouse},prt_conditions,mouse{i_mouse})
        %         end
    end
end

%% Plot mean anticipatory licking in all mice and mean across mice

[~, ind_CS, ~, ~] = get_indicesEvents(time_licks,fluotype,CSin,1,2.8,TraceEnd);
[~, ind_CS_rotenc, ~, ~] = get_indicesEvents(time_rotenc,fluotype,CSin,1,2.8,TraceEnd);

if num_CSTT > 0 && taskPhase ~= 2 % do not do this plot in the case of reversal learning
    alllicks = nan(numel(time_licks),num_mice,num_CSTT);
    for i_TT = 1:num_CSTT
        for i_mouse = 1:num_mice
            alllicks(:,i_mouse,i_TT) = lick_mean_timecourse{i_mouse}{i_TT};
        end
    end
    alllicks_means = squeeze(mean(alllicks(ind_CS,:,:)));

    % plot mean lick trace along time
    figure('Position',[50 100 500 300]); hold on
    fill([0 1 1 0],[0 0 10 10],[0.6 0.6 0.6],'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
    line([3 3],[0 10],'color','k')
    for i_TT = 1:num_CSTT
        erro = nanstd(alllicks(:,:,i_TT),[],2)/sqrt(size(alllicks,2));
        media = nanmean(alllicks(:,:,i_TT),2);
        pos1 = find(~isnan(erro),1,'first');
        pos2 = find(~isnan(erro),1,'last');
        fill([time_licks(pos1:pos2) fliplr(time_licks(pos1:pos2))],[media(pos1:pos2)+erro(pos1:pos2);flipud(media(pos1:pos2)-erro(pos1:pos2))]',prt_colors(i_TT,:),'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
        plot(time_licks,media,'color',prt_colors(i_TT,:),'linewidth',2)
    end
    set(gca,'Tickdir','out','TickLength',2*(get(gca,'TickLength')),'Fontsize',12)
    xlim([-1 6])
    ylim([0 7])
    xlabel('Time (s)')
    ylabel('Licks/s')

    % plot mean anticipatoy licking for each trial type as inset
    axes('Position',[.7 .7 .2 .2]); hold on
    box on
    plot(1:num_CSTT,alllicks_means,'color',[.6 .6 .6])
    for i = 1:num_CSTT
        erro = std(alllicks_means(:,i)/sqrt(num_mice));
        errorbar(i,mean(alllicks_means(:,i)),erro,'color',prt_colors(i,:),'linewidth',2) 
        plot(i,mean(alllicks_means(:,i)),'o','MarkerFaceColor',prt_colors(i,:),'MarkerEdgeColor',prt_colors(i,:),'MarkerSize',5)
    end
    set(gca,'Tickdir','out','TickLength',2*(get(gca,'TickLength')),'xlim',[0.5 num_CSTT+0.5])
    xlabel('Trial type')
%     ylabel('Licks/s')
    title('Mean anticipatory licking')
    xticks(1:6)
    xticklabels(prt_conditions(1:num_CSTT))
    
    % For rotary encoder
    allrotenc = nan(numel(time_rotenc),num_mice,num_CSTT,size(rotencvar,2));
    for i_TT = 1:num_CSTT
        for i_mouse = 1:num_mice
            allrotenc(:,i_mouse,i_TT,:) = rotenc_mean_timecourse{i_mouse}{i_TT};
        end
    end
    allrotenc_means = squeeze(mean(allrotenc(ind_CS_rotenc,:,:,:)));
    
    % plot mean rotary encoder traces along time
    for i_rotencvar = 1:size(rotencvar,2)        
        figure('Position',[50 100 500 300]); hold on
        fill([0 1 1 0],[-10 -10 10 10],[0.6 0.6 0.6],'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
        line([3 3],[0-1 10],'color','k')
        for i_TT = 1:num_CSTT
            erro = nanstd(allrotenc(:,:,i_TT,i_rotencvar),[],2)/sqrt(size(allrotenc,2));
            media = nanmean(allrotenc(:,:,i_TT,i_rotencvar),2);
            pos1 = find(~isnan(erro),1,'first');
            pos2 = find(~isnan(erro),1,'last');
            fill([time_rotenc(pos1:pos2) fliplr(time_rotenc(pos1:pos2))],[media(pos1:pos2)+erro(pos1:pos2);flipud(media(pos1:pos2)-erro(pos1:pos2))]',prt_colors(i_TT,:),'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
            plot(time_rotenc,media,'color',prt_colors(i_TT,:),'linewidth',2)
        end
        set(gca,'Tickdir','out','TickLength',2*(get(gca,'TickLength')),'Fontsize',12)
        xlim([-1 6])
        xlabel('Time (s)')
        ylabel(rotencvar(i_rotencvar))
                ylim([-0.5 0.5])
    
        % plot mean anticipatoy rotary encoder info for each trial type as inset
        axes('Position',[.7 .7 .2 .2]); hold on
        box on
        plot(1:num_CSTT,allrotenc_means(:,:,i_rotencvar),'color',[.6 .6 .6])
        for i = 1:num_CSTT
            erro = std(allrotenc_means(:,i,i_rotencvar)/sqrt(num_mice));
            errorbar(i,mean(allrotenc_means(:,i,i_rotencvar)),erro,'color',prt_colors(i,:),'linewidth',2) 
            plot(i,mean(allrotenc_means(:,i,i_rotencvar)),'o','MarkerFaceColor',prt_colors(i,:),'MarkerEdgeColor',prt_colors(i,:),'MarkerSize',5)
        end
        set(gca,'Tickdir','out','TickLength',2*(get(gca,'TickLength')),'xlim',[0.5 num_CSTT+0.5])
        xlabel('Trial type')
    %     ylabel('Licks/s')
        title(['Mean anticipatory ' rotencvar(i_rotencvar)])
        xticks(1:6)
        xticklabels(prt_conditions(1:num_CSTT))
    end
end 

%% Formal Analysis of Neural Data

if taskPhase == 1 || taskPhase == 3
    % Distributional RL analsys
    run_analysis_DRL_MSFP(protocol,fluotype,mouse,day,num_TT,num_cells,fluo_all,time_fluo,RotEnc_all,CSin,TraceEnd,revtrial,SR_fluo,Rsize_pos,TT_allRsizes,fiberIDs)
elseif taskPhase == 0 || taskPhase == 2
    % Learning analysis: either initial learning or reversal learning
    run_analysis_Learning_MSFP(protocol,taskPhase,TT_shift,fluotype,smooth_acrosstrials,mouse,day,num_TT,num_cells,licks_all,time_licks,fluo_all,time_fluo,ind_CS,CSin,TraceEnd,revtrial,fiberIDs)
end

end





