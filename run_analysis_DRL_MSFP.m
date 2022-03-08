function run_analysis_DRL_MSFP(protocol,fluotype,mouse,day,num_TT,num_cells,fluo_all,time_fluo,rotenc_all,CSin,TraceEnd,revtrial,SR_fluo,Rsize_pos,TT_allRsizes,fiberIDs)

[prt_colors, prt_conditions, num_CSTT, CSTT_mostR, R_colors] = get_protocol_info(protocol);
num_mice = length(mouse);


%% Perform an analysis to compare different regions in terms of CS and US responses
num_ROIs = num_cells;
mat_analysis_Regions_MSFP(protocol,fluotype,mouse,day,num_TT,num_ROIs,fluo_all,time_fluo,rotenc_all,CSin,TraceEnd,SR_fluo,Rsize_pos,TT_allRsizes,fiberIDs);

%% Separate trials of the same type with different rewards, Re-organize cells in matrices
% Mean US response to all trials (dimensions: cell id, reward size id,
% trials id (put nans when the number of trials is not always the same for
% all conditions)
% create matrices to feed in zero point crossing calculation from DeepMind
% fluo_all_org is a cell with an entry for each TT, the dim of each cell are: time,cellid,Rsize,trial#
% fluo_all_org_info

fluo_all_org_temp = cell(1,max(num_TT));
cellid = cell(1,sum(num_cells));
ind = 1;

for i_TT = 1:max(num_TT)
    
    max_numTrials = nan(1,num_mice);
    num_Rsize = max_numTrials;
    for i_mouse = 1:num_mice
        if i_TT <= max(num_TT(i_mouse)) % some mice do not have TT7
            % max number of trials for reward sizes delivered in this TT
            max_numTrials(i_mouse) = max(cellfun(@(x) size(x,1),Rsize_pos{i_mouse}(i_TT,:),'uniform',true));
            % max number of R sizes available in this TT
            num_Rsize(i_mouse) = length(TT_allRsizes{i_mouse}{i_TT});
        end
    end
    max_numTrials = max(max_numTrials);
    max_numRsizes = max(num_Rsize);
    
    for i_mouse = 1:num_mice
        if i_TT <= max(num_TT(i_mouse))
            if i_mouse == 1  % allocate matrix (cells of all mice will be concatenated)
                fluo_all_org_temp{i_TT} = nan(length(time_fluo),sum(num_cells),max_numRsizes,max_numTrials);
            end
            for i_Rsize = 1:num_Rsize(i_mouse)
                % trial positions for this R size of this TT in this mouse
                pos_Rsize =  Rsize_pos{i_mouse}{i_TT,i_Rsize};
                data =  fluo_all{i_mouse}{i_TT}(:,pos_Rsize,:);
                
                num_trialsRsize = size(data,2);
                if i_mouse == 1
                    indx1_cells = 1;
                else
                    indx1_cells = sum(num_cells(1:i_mouse-1))+1;
                end
                indx2_cells = sum(num_cells(1:i_mouse));
                
                fluo_all_org_temp{i_TT}(:,indx1_cells:indx2_cells,i_Rsize,1:num_trialsRsize) = permute(data,[1 3 4 2]);
                %                 display(['num 0s in data TT' num2str(i_TT) 'Rsize' num2str(i_Rsize) ':' num2str(length(find(data(:)==0)))])
                %                 display(['num 0s in fluoorg:' num2str(length(find(fluo_all_org_temp{i_TT}(:)==0)))])
            end
        end
        if i_TT == 1
            % save info about cellid
            for i_cell = 1:num_cells(i_mouse)
                cellid{ind} = [mouse{i_mouse} '_' day{i_mouse}, '_' fiberIDs{i_cell}];
                ind = ind+1;
            end
            % save Rsizes available
            if i_mouse == num_mice  % this is an examples mouse with all possible rewards
                Rsizes =  TT_allRsizes{i_mouse};
            end
            %             if strcmp(mouse(i_mouse),'SM104')  % this is an examples mouse with all possible rewards
            %                 Rsizes =  TT_allRsizes{i_mouse};
            %             end
        end
    end
end

% Not entirely clear: when using permute, sometimes 0s are added to the
% matrices, so I will change those to nans (because they will bias my
% averages later): this actually is happening only in TT7, probably because
% some mice do not have trials there, but not sure why this happens
fluo_all_org = fluo_all_org_temp;
for i_TT = 1:max(num_TT)
    fluo_all_org{i_TT}(fluo_all_org_temp{i_TT} == 0) = nan;
end

%% Plot mean of all cells and check that we get what is expected

plot_pcolorplot_meanallcells(fluo_all_org,time_fluo,prt_conditions,...
    num_CSTT,CSin{1},TraceEnd{1},revtrial{1},SR_fluo{1},1,Rsizes);

%% Define periods for CS and US

CSstartT = 0.15;  % quantify from 100 ms after CS start
CSendT = 1.5;  % quantify until 500 ms before US start
[ind_BSL, ind_CS, ind_US_exp, ind_US_unexp, ind_US_exp_early] = get_indicesEvents(time_fluo,fluotype,CSin,CSstartT,CSendT,TraceEnd);

baselines = cellfun(@(x) nanmean(x(ind_BSL,:,:,:),1),fluo_all_org,'UniformOutput',false);  % should not squeeze, because I don't want to squeeze the 3rd dimension
CStrans = cellfun(@(x) nanmean(x(ind_CS,:,:,:),1),fluo_all_org,'UniformOutput',false);
UStrans = cellfun(@(x) nanmean(x(ind_US_exp,:,:,:),1),fluo_all_org,'UniformOutput',false);
if num_TT > num_CSTT % if there are trials with unexpected reward delivery
    UStrans{num_CSTT+1} = nanmean(fluo_all_org{num_CSTT+1}(ind_US_unexp,:,:,:),1);
    CStrans{num_CSTT+1} = nan(size(CStrans{num_CSTT+1}));
end

% To compare CS responses to baseline, concatenate all possible R sizes in the TT
baselines_compCS = cellfun(@(x) reshape(x,size(x,1),size(x,2),size(x,3)*size(x,4)),baselines,'UniformOutput',false);
CStrans_compBSL = cellfun(@(x) reshape(x,size(x,1),size(x,2),size(x,3)*size(x,4)),CStrans,'UniformOutput',false);

%% Formal DRL Analysis:
roi_selection_method = 'None';
switch roi_selection_method
    case 'None'
        fluo_all_org_sel = fluo_all_org;
        baselines_sel = baselines;
        CStrans_sel = CStrans;
        UStrans_sel = UStrans;
        cellid_sel = cellid;
end
%% Compare US transients to expected and unexpected rewards

mat_analysis_US(protocol,cellid_sel,UStrans_sel,Rsizes,fluotype,...
    prt_conditions, num_CSTT, CSTT_mostR, R_colors)

%% Compare CS responses, their variance and classify neurons

mat_analysis_CS(cellid_sel,CStrans_sel,Rsizes,fluotype,protocol)

%% Determine Reversal Points and compare two halves of the data

ZCmethod = 'sum'; % options: 'count', *'sum'*, 'proportion', ('activity' this last one is not working)
ZCpoints_all = mat_analysis_ZC(cellid_sel,UStrans_sel,Rsizes,fluotype,protocol,ZCmethod);

%% Calculate asymmetries

functionType = 'linear'; % options: 'linear', 'hill'
mat_analysis_Asym(cellid_sel,UStrans_sel,Rsizes,fluotype,protocol,...
    ZCmethod,functionType,CStrans_sel)

%% Relationship between zero-crossing and asymmetric scaling

[TausH, TausT, ZCpointsH] = mat_analysis_AsymVsZC(cellid_sel,UStrans_sel,Rsizes,fluotype,protocol,...
    ZCmethod,functionType);
% correr esta funcao 3 vezes para ter nice plots

switch protocol
    case 'DRL6Od'
        sZCpoints_TT2 = sprintf('%.15f,',ZCpointsH{1});
        sZCpoints_TT2 = sZCpoints_TT2(1:end-1);
        sTaus_TT2 = sprintf('%.15f,',TausH{1});
        sTaus_TT2 = sTaus_TT2(1:end-1);
        
        sZCpoints_TT3 = sprintf('%.15f,',ZCpointsH{2});
        sZCpoints_TT3 = sZCpoints_TT3(1:end-1);
        sTaus_TT3 = sprintf('%.15f,',TausH{2});
        sTaus_TT3 = sTaus_TT3(1:end-1);
        
        sZCpoints_TT4 = sprintf('%.15f,',ZCpointsH{3});
        sZCpoints_TT4 = sZCpoints_TT4(1:end-1);
        sTaus_TT4 = sprintf('%.15f,',TausH{3});
        sTaus_TT4 = sTaus_TT4(1:end-1);
        
    case 'VarMag'
        sZCpoints = sprintf('%.15f,',ZCpointsH{1});
        sZCpoints = sZCpoints(1:end-1);
        
        sTaus = sprintf('%.15f,',TausH{1});
        sTaus = sTaus(1:end-1);
        
        sTaus_T = sprintf('%.15f,',TausT{1});
        sTaus_T = sTaus_T(1:end-1);
        
    case 'DRL6Rv'
        sZCpoints_TT5 = sprintf('%.15f,',ZCpointsH{1});
        sZCpoints_TT5 = sZCpoints_TT5(1:end-1);
        sTaus_TT5 = sprintf('%.15f,',TausH{1});
        sTaus_TT5 = sTaus_TT5(1:end-1);
        
        sZCpoints_TT6 = sprintf('%.15f,',ZCpointsH{2});
        sZCpoints_TT6 = sZCpoints_TT6(1:end-1);
        sTaus_TT6 = sprintf('%.15f,',TausH{2});
        sTaus_TT6 = sTaus_TT6(1:end-1);
end


%% Location analysis

if strcmp(protocol,'DRL6Od') || strcmp(protocol,'VarMag') || strcmp(protocol,'Eshl16') %% This need to be generalized to other protocols
    ASYMmethod = {'alpha', 'theta'}; % options: 'theta', 'alpha'; Was not defined before because the above functions calculate and plot both
    ASYMmethod = {'alpha'};
    for i = 1:length(ASYMmethod)
        [Taus{i}, ZCpoints{i}, CellIDTaus{i},ZCpoints_TT4] = mat_analysis_OptimismLocation_MSFP(cellid_sel,UStrans_sel,Rsizes,fluotype,...
            protocol,ZCmethod,ASYMmethod{i},functionType);
    end
end

% sTaus = sprintf('%.15f,',Taus{1});
% sTaus = sTaus(1:end-1);
% sZCpoints = sprintf('%.15f,',ZCpoints{1});
% sZCpoints = sZCpoints(1:end-1);
%
% sZCpoints
% sTaus

% if strcmp(protocol,'VarMag')
%     savefolder = '/home/sara/Documents/DATA/Imaging/Analysis/Aug2020/';
%     savename = ['VarMag_Taus_' fluotype '_' date];
%     save([savefolder savename],'Taus','ZCpoints','CellIDTaus')
% end

%% Relationship between CS and US responses, with info acount asymmtries and ZC

% This only works in Sara's Uchida Lab Desktop because it's where I
% have the files with masks
mat_analysis_CSUS(protocol, cellid_sel,CStrans_sel,UStrans_sel, Rsizes,...
    fluotype,R_colors,Taus, ZCpoints, CellIDTaus,ZCpoints_TT4)

%% Organize asymmetries and reversal points for each neuron: prepare for decoding

if strcmp(protocol,'DRL6Od') || strcmp(protocol,'VarMag')
    % ZCpoints_all: ZC of all active cells in this task, for TT2,3,4
    % cellid_act: ids of all active cells (match the ZCpoints_all values)
    % Taus: taus of active cells that could be fit for asymmetry
    % calculation. Includes both methods of calculating taus (alpha and
    % theta)
    % ZCpoints: ZC points in TT4 of cells from Taus
    % CellIDTaus: cell id of cells from Taus
    prepareDecoding(protocol,fluotype,ZCpoints_all, cellid_sel, Taus, ZCpoints, CellIDTaus)
end

%% PCA analysis

%% Tensor analysis

%% GLM analysis

end