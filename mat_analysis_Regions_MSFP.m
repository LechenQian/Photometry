function mat_analysis_Regions_MSFP(protocol,fluotype,mouse,day,num_TT,num_ROIs,fluo_all,time_fluo,rotenc_all,CSin,TraceEnd,SR_fluo,Rsize_pos,TT_allRsizes,fiberIDs)

[prt_colors, prt_conditions, num_CSTT, CSTT_mostR, R_colors] = get_protocol_info(protocol);

%% This is needed for the 'fill' function 
pos1 = find(time_fluo >= -1.8, 1, 'first');
pos2 = find(time_fluo <= 8, 1,'last');
time_toplot = time_fluo(pos1:pos2);

num_mice = numel(mouse);
num_TT = num_TT(1);
num_time = numel(time_fluo);
num_ROI = num_ROIs(1);
num_Rsz = numel(TT_allRsizes{1}{CSTT_mostR});

%% Calculate the z-scored fluorescence: z-scoring each ROI of each mouse across the entire session

fluo_all_zscore = cell(1,num_mice);

for i_mouse = 1:num_mice
    for i_ROI = 1:num_ROI
        
        vectorSizes1 = cellfun(@(x) [size(x,1) size(x,2)],fluo_all{i_mouse},'UniformOutput',false);
        
        fluo_forzscore = cellfun(@(x) reshape(squeeze(x(:,:,i_ROI)),numel(x(:,:,i_ROI)),1),fluo_all{i_mouse},'UniformOutput',false);
        vectorSizes2 = cellfun(@(x) numel(x),fluo_forzscore);
        
        fluo_forzscore_mat = cell2mat(fluo_forzscore');
        fluo_zscore_aux = normalize(fluo_forzscore_mat);
        
        fluo_zscore_cell = mat2cell(fluo_zscore_aux,vectorSizes2,1)';
        
        for i_TT = 1:num_TT
            fluo_all_zscore{i_mouse}{i_TT}(:,:,i_ROI) = reshape(fluo_zscore_cell{i_TT},vectorSizes1{i_TT});
        end
    end
end

%% Plot fluorescence timecourse across different regions (deltaF/F and z-score)

fluos = {'fluo_all','fluo_all_zscore'};
fluos_labels = {'\DeltaF/F(%)','z-score'};
fluomult = [100 1];

for i_fluos = 1:numel(fluos)
    fluo_toplot = eval(fluos{i_fluos});
    fluo_mean_timecourse = cell(1,num_mice);
    fluo_std_timecourse =  cell(1,num_mice);
    for i_mouse = 1:num_mice
        fluo_mean_timecourse{i_mouse} = cellfun(@(x) squeeze(nanmean(x,2)),fluo_toplot{i_mouse},'UniformOutput',false);
        fluo_std_timecourse{i_mouse} = cellfun(@(x) squeeze(nanstd(x,0,2)),fluo_toplot{i_mouse}, 'UniformOutput',false);
    end

    % Make one figure per trialtype
    for i_TT = 1:num_TT
        figure('Position', [200 50 1000 500],'name',[fluos_labels{i_fluos} ' TT' num2str(i_TT)],'numbertitle','off', 'Units', 'Pixels');
        set(gcf , 'Color', 'w'); hold on
        for i_ROI = 1:num_ROI
           h = subplot(3,num_ROI/3,i_ROI); hold on
           for i_mouse = 1:num_mice   
               paramedia(:,i_mouse) = fluomult(i_fluos)*fluo_mean_timecourse{i_mouse}{i_TT}(:,i_ROI);
               plot(time_fluo,paramedia(:,i_mouse),'color',[1 1 1]*0.2*i_mouse)            
           end
           media = nanmean(paramedia(pos1:pos2,:),2);
           erroaux = nanstd(paramedia(pos1:pos2,:),[],2)/sqrt(num_mice);
           mediaerro = media(~isnan(erroaux));
           erro = erroaux(~isnan(erroaux));
           time_toplot = time_toplot(~isnan(erroaux));               

           fill([time_toplot fliplr(time_toplot)],[mediaerro+erro; flipud(mediaerro-erro)],[1 0.2 0.5],'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
           plot(time_toplot,mediaerro,'color',[1 0.2 0.5],'linewidth',2) 
           title(fiberIDs{i_ROI})
           ylabel(fluos_labels{i_fluos})
           if i_fluos == 1
            ylim([-2 4])
           else
            ylim([-0.8 3])
           end
           xlabel('Time - Cue (s)')
           xlim([-2 6])
           set(h,'TickDir','out')
        end
    end
end

%% Do the same plot with averages separated by R size

mean_Rsize_permouse = nan(num_TT,num_mice,num_Rsz,num_time,num_ROI);
mean_Rsize = nan(num_TT,num_Rsz,num_time,num_ROI);
mean_Rsize_sem = mean_Rsize;

for i_TT = 1:num_TT
    figure('Position', [200 50 1000 500],'name',[fluos_labels{i_fluos} ' TT' num2str(i_TT)],'numbertitle','off', 'Units', 'Pixels');
    set(gcf , 'Color', 'w'); hold on
    num_Rsz_TT = numel(TT_allRsizes{1}{i_TT});
    for i_rsz = 1:num_Rsz_TT
        for i_mouse = 1:num_mice
            pos_thisRsz  = Rsize_pos{i_mouse}{i_TT,i_rsz};
            mean_Rsize_permouse(i_TT,i_mouse,i_rsz,:,:) = squeeze(nanmean(fluo_all_zscore{i_mouse}{i_TT}(:,pos_thisRsz,:),2));
        end
        mean_Rsize(i_TT,i_rsz,:,:) = squeeze(nanmean(mean_Rsize_permouse(i_TT,:,i_rsz,:,:),2));
        mean_Rsize_sem(i_TT,i_rsz,:,:) = squeeze(nanstd(mean_Rsize_permouse(i_TT,:,i_rsz,:,:),[],2)/sqrt(num_mice));
        rsz_color = R_colors{i_TT}(i_rsz,:);
        
        for i_ROI = 1:num_ROI
            h = subplot(3,num_ROI/3,i_ROI); hold on
            media = squeeze(mean_Rsize(i_TT,i_rsz,pos1:pos2,i_ROI));
            erroaux = squeeze(mean_Rsize_sem(i_TT,i_rsz,pos1:pos2,i_ROI));
            
            mediaerro = media(~isnan(erroaux));
            erro = erroaux(~isnan(erroaux));
            time_toplot = time_toplot(~isnan(erroaux));
            
%             fill([time_toplot fliplr(time_toplot)],[mediaerro+erro; flipud(mediaerro-erro)],rsz_color,'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
            plot(time_toplot,mediaerro,'color',rsz_color,'linewidth',1)
            title(fiberIDs{i_ROI})
            ylabel('z-score')
            ylim([-0.8 4])
            xlabel('Time - Cue (s)')
            xlim([-2 6])
            set(h,'TickDir','out')
        end
    end
end 
        

%% Do one subplot per reward size and plot all regions together

% Set color for each region:
% FiberColors = tab20(num_ROI-1);
% FiberColors = brewermap(num_ROI-1,'*PiYG');
FiberColors = brewermap(num_ROI-1,'*Spectral');


for i_TT = 1:num_TT
    figure('Position', [200 50 1000 500],'name',[fluos_labels{i_fluos} ' TT' num2str(i_TT)],'numbertitle','off', 'Units', 'Pixels');
    set(gcf , 'Color', 'w'); hold on
    num_Rsz_TT = numel(TT_allRsizes{1}{i_TT});
    for i_rsz = 1:num_Rsz_TT
        h = subplot(3,ceil(num_Rsz/3),i_rsz); hold on
        
        media = squeeze(mean_Rsize(i_TT,i_rsz,pos1:pos2,1:end-1));
        
        ho = plot(time_toplot,media);
        title([num2str(TT_allRsizes{1}{i_TT}(i_rsz)) '\mul'])
        ylabel('z-score')
        ylim([-0.8 4])
        xlabel('Time - Cue (s)')
        xlim([-2 6])
        set(h,'TickDir','out')
        set(ho,{'color'},num2cell(FiberColors,2))
    end
    subplot(3,ceil(num_Rsz/3),i_rsz+1);
    h2 = plot(time_toplot,media);
    set(h2,{'color'},num2cell(FiberColors,2))
    legend(fiberIDs,'NumColumns',3)
    set(gca,'visible','off')
end


%% Plot mean transient to different USs for each region

% Calculate CS and US transients

CSstartT = 0.1;  % quantify from 100 ms after CS start
CSendT = 1.0;  % quantify until 500 ms before US start
[ind_BSL, ind_CS, ind_US_exp, ind_US_unexp, ind_US_exp_early] = get_indicesEvents(time_fluo,fluotype,CSin,CSstartT,CSendT,TraceEnd);

% mean_Rsize has dim: TT, Rsz, time, ROI
% mean_Rsize_permouse dim: TT, mouse,rsz,time,ROI

baselines = squeeze(nanmean(mean_Rsize_permouse(:,:,:,ind_BSL,:),4));
CStrans = squeeze(nanmean(mean_Rsize_permouse(:,:,:,ind_CS,:),4));
UStrans_early = squeeze(nanmean(mean_Rsize_permouse(:,:,:,ind_US_exp_early,:),4));
UStrans = squeeze(nanmean(mean_Rsize_permouse(:,:,:,ind_US_exp,:),4));
if num_TT > num_CSTT % if there are trials with unexpected reward delivery
    UStrans(num_CSTT+1,:,:,:) = squeeze(nanmean(mean_Rsize_permouse(num_CSTT+1,:,:,ind_US_unexp,:),4));
    CStrans(num_CSTT+1,:,:,:) = nan(size(CStrans(num_CSTT+1,:,:,:)));
end
num_Rsz_max = size(UStrans,3);

% To compare CS responses to baseline, concatenate all possible R sizes in the TT
% baselines_compCS = cellfun(@(x) reshape(x,size(x,1),size(x,2),size(x,3)*size(x,4)),baselines,'UniformOutput',false);
% CStrans_compBSL = cellfun(@(x) reshape(x,size(x,1),size(x,2),size(x,3)*size(x,4)),CStrans,'UniformOutput',false);

% Plot US amplitudes:
for i_TT = 1:num_TT
    fUS(i_TT) = figure('Position', [200 50 1000 500],'name',['US amplitudes - TT' num2str(i_TT)],'numbertitle','off', 'Units', 'Pixels');
    set(gcf , 'Color', 'w'); hold on
    for i_ROI = 1:num_ROI
       h = subplot(3,ceil(num_ROI/3),i_ROI); hold on
       plot([5 5],[-0.8 3],'--','color','k','HandleVisibility','off')
       plot([0 max(TT_allRsizes{1}{CSTT_mostR})+1],[0 0],'--','color','k','HandleVisibility','off')
       plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(UStrans_early(i_TT,:,:,i_ROI)),':','color',[0.6 0.6 0.6],'HandleVisibility','off');
       plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(mean(UStrans_early(i_TT,:,:,i_ROI),2)),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2) ;
       
       plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(UStrans(i_TT,:,:,i_ROI)),'color',[0.6 0.6 0.6],'HandleVisibility','off') ;
       plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(mean(UStrans(i_TT,:,:,i_ROI),2)),'-o','MarkerFaceColor',[1 0.2 0.5],'color',[1 0.2 0.5],'linewidth',2) ;
       
       title(fiberIDs{i_ROI})
       if i_ROI == 1
           ylabel('z-score')
       end
       ylim([-0.8 3])
       if i_ROI == num_ROI
           xlabel('Reward (\mul)','Interpreter','latex');
           legend('Early','Late')
       end
       xticks(TT_allRsizes{1}{CSTT_mostR}');
       set(gca,'TickDir','out','xticklabels',cellstr(num2str(TT_allRsizes{1}{i_TT}))')
    end
end

%% Fitting curves to get the reponse function and ZC for each region:

% Using  the same function that Mitsuko used in Iku's first eLife paper
% does not work well in this data set because it does not assimptote
% model_fun = fittype('a*(x^0.7+c)'); 

% model_fun_hill = fittype('beta + fmax*((x.^0.5)./(sigma.^0.5 + x.^0.5))');
% model_fun_linear = fittype('beta + gamma*x');

model_fun_hill = fittype(@(beta,fmax,sigma,x) beta + fmax*((x.^0.5) ./ (sigma.^0.5 + x.^0.5)));
model_fun_linear = fittype(@(beta,gamma,x) beta + gamma*x);

xfit = TT_allRsizes{1}{CSTT_mostR};
UStrans_fitted = nan(num_mice,numel(xfit),num_ROI);
UStrans_fitted_early = nan(num_mice,numel(xfit),num_ROI);

figure('Position', [200 50 1000 500],'name',['US amplitudes - TT' num2str(CSTT_mostR)],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
for i_ROI = 1:num_ROI             % find fitting parameters for each region
    for i_mouse = 1:num_mice
        yfit = squeeze(UStrans(CSTT_mostR,i_mouse,:,i_ROI));
        try
            [fitobj{i_mouse,i_ROI},gof{i_mouse,i_ROI},...
                output{i_mouse,i_ROI}] = fit(xfit,yfit_early,model_fun_hill,'StartPoint',[0; 50; 0.5]);
            beta = fitobj{i_mouse,i_ROI}.beta;
            fmax = fitobj{i_mouse,i_ROI}.fmax;
            sigma = fitobj{i_mouse,i_ROI}.sigma;
            UStrans_fitted(i_mouse,:,i_ROI) = model_fun_hill(beta,fmax,sigma,xfit);
            
            func = @(x) beta + fmax*((x.^0.5)./(sigma.^0.5 + x.^0.5));
            x0 = 1;
            x0_calc = fzero(func,x0);
            ZC(i_mouse,i_ROI) = x0_calc;
        catch
            try
                [fitobj{i_mouse,i_ROI},gof{i_mouse,i_ROI},...
                    output{i_mouse,i_ROI}] = fit(xfit,yfit_early,model_fun_linear,'StartPoint',[0; 1]);
                beta = fitobj{i_mouse,i_ROI}.beta;
                gamma = fitobj{i_mouse,i_ROI}.gamma;
                UStrans_fitted(i_mouse,:,i_ROI) = model_fun_linear(beta,gamma,xfit);
                
                func = @(x) beta_early + gamma_early*x;
                x0 = 1;
                x0_calc = fzero(func,x0);
                ZC(i_mouse,i_ROI) = x0_calc;
            catch
                UStrans_fitted_early(i_mouse,:,i_ROI) = nan(size(UStrans_fitted(i_mouse,:,i_ROI)));
                ZC(i_mouse,i_ROI) = nan;
            end
        end
        
        yfit_early = squeeze(UStrans_early(CSTT_mostR,i_mouse,:,i_ROI));
        try            
            [fitobj_early{i_mouse,i_ROI},gof_early{i_mouse,i_ROI},...
                output_early{i_mouse,i_ROI}] = fit(xfit,yfit_early,model_fun_hill,'StartPoint',[0; 50; 0.5]);
            beta_early = fitobj_early{i_mouse,i_ROI}.beta;
            fmax_early = fitobj_early{i_mouse,i_ROI}.fmax;
            sigma_early = fitobj_early{i_mouse,i_ROI}.sigma;
            UStrans_fitted_early(i_mouse,:,i_ROI) = model_fun_hill(beta_early,fmax_early,sigma_early,xfit);
            
            func =@(x) beta_early + fmax_early*((x.^0.5)./(sigma_early.^0.5 + x.^0.5));
            x0 = 1;
            x0_calc = fzero(func,x0);
            ZC_early(i_mouse,i_ROI) = x0_calc;
        catch
%             try
%                [fitobj_early{i_mouse,i_ROI},gof_early{i_mouse,i_ROI},...
%                 output_early{i_mouse,i_ROI}] = fit(xfit,yfit_early,model_fun_linear,'StartPoint',[0; 1]);
%                 beta_early = fitobj_early{i_mouse,i_ROI}.beta;
%                 gamma_early = fitobj_early{i_mouse,i_ROI}.gamma;
%                 UStrans_fitted_early(i_mouse,:,i_ROI) = model_fun_linear(beta_early,gamma_early,xfit); 
%                 
%                 func = @(x) beta_early + gamma_early*x;                
%                 x0 = 1;
%                 x0_calc = fzero(func,x0);
%                 ZC_early(i_mouse,i_ROI) = x0_calc;
%             catch
                UStrans_fitted_early(i_mouse,:,i_ROI) = nan(size(UStrans_fitted_early(i_mouse,:,i_ROI))); 
                ZC_early(i_mouse,i_ROI) = nan;
%             end
        end
    end
    
    h = subplot(3,ceil(num_ROI/3),i_ROI); hold on
    plot([5 5],[-0.8 3],'--','color','k','HandleVisibility','off')
    plot([0 max(TT_allRsizes{1}{CSTT_mostR})+1],[0 0],'--','color','k','HandleVisibility','off')
    
    plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(UStrans_fitted_early(:,:,i_ROI)),':','color',[0.6 0.6 0.6],'HandleVisibility','off');
    plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(nanmean(UStrans_fitted_early(:,:,i_ROI))),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2) ;
    text(1,2,['ZCerl=' num2str(nanmean(ZC_early(:,i_ROI)))]);
    
    plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(UStrans_fitted(:,:,i_ROI)),':','color',[0.6 0.6 0.6],'HandleVisibility','off');
    plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(nanmean(UStrans_fitted(:,:,i_ROI))),'-o','MarkerFaceColor',[1 0.2 0.5],'color',[1 0.2 0.5],'linewidth',2) ;
    text(1,1,['ZC=' num2str(nanmean(ZC(:,i_ROI)))]);
    
    title(fiberIDs{i_ROI})
    if i_ROI == 1
        ylabel('z-score')
    end
    ylim([-0.8 3])
    if i_ROI == num_ROI
        xlabel('Reward (\mul)','Interpreter','latex');
        legend('Early','Late')
    end
    xticks(TT_allRsizes{1}{CSTT_mostR}');
    set(gca,'TickDir','out','xticklabels',cellstr(num2str(TT_allRsizes{1}{i_TT}))')
end

figure(fUS(CSTT_mostR))
for i_ROI = 1:num_ROI  
    subplot(3,ceil(num_ROI/3),i_ROI); hold on
    plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(nanmean(UStrans_fitted_early(:,:,i_ROI))),'color',[1 0.5 0.8]-0.5,'linewidth',2) ;
    plot(TT_allRsizes{1}{CSTT_mostR}',squeeze(nanmean(UStrans_fitted(:,:,i_ROI))),'color',[1 0.2 0.5]-0.2,'linewidth',2) ;
end
 
% Plot zero-crossings
figure('Position', [200 50 500 500],'name','Zero-crossing','numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
plot(1:num_ROI,ZC_early,':','color',[0.6 0.6 0.6],'HandleVisibility','off') ;
plot(1:num_ROI,nanmean(ZC_early),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2);
plot(1:num_ROI,ZC,'color',[0.6 0.6 0.6],'HandleVisibility','off') ;
plot(1:num_ROI,nanmean(ZC),'-o','MarkerFaceColor',[1 0.2 0.5],'color',[1 0.2 0.5],'linewidth',2);
ylabel('ZC (from fit)')
legend('Early','Late')
xticks(1:num_ROI);
xtickangle(45)
set(gca,'TickDir','out','xticklabels',fiberIDs)

%% Plot these zero-crossings in 3D map

% Coordinates of the bottom of the leftmost shank (channels facing us)
MCoor = table_MouseImplantFibers_MSFP;
for i_cell = 1:num_ROIs
    ROI = fiberIDs{i_cell}; 
    AP(i_cell) = table2array(MCoor(ROI,'AP_target'));              % (1:end-4) is just to remove the ' .mat '
    ML(i_cell) = table2array(MCoor(ROI,'ML_target'));
    DV(i_cell) = table2array(MCoor(ROI,'DV_target'));
%     hemisph = char(table2array(ROI,'ImplantHemisphere'));
end

figure('Position',[50 50 750 400]);
LocAxis = cat(1,AP,ML,DV);
LocAxislabel = {'AP (mm)','ML (mm)', 'DV (mm)'};
Locylim = [-2 2; 0 4; -5 -1];

    
for i = 1:size(LocAxis,1)
    subplot(2,3,i); hold on
    for i_cell = 1:num_ROIs-4
        plot(LocAxis(i,i_cell),ZC(:,i_cell),'o')
    end 
    LocAxis_reg = repmat(LocAxis(i,1:end-4),num_mice,1);
    ZC_reg = ZC(:,1:end-4);
    ZC_reg = ZC_reg(:);
    ZC_reg = ZC_reg(~isnan(ZC_reg));
    LocAxis_reg = LocAxis_reg(:);
    LocAxis_reg = LocAxis_reg(~isnan(ZC_reg));
    LocAxis_reg = LocAxis_reg(~isnan(LocAxis_reg));
    ZC_reg = ZC_reg(~isnan(LocAxis_reg));
    
    [beta, ~, stats] = glmfit(LocAxis_reg,ZC_reg);
    R = corrcoef(LocAxis_reg,ZC_reg);
    plot([min(LocAxis_reg)-0.1 max(LocAxis_reg)+0.1], ...
        beta(1) + beta(2)*[min(LocAxis_reg)-0.1 max(LocAxis_reg)+0.1], 'k-')
    title(['R = ' num2str(round(R(1,2),2)) ', p = ' num2str(round(stats.p(2),4))])
    xlabel(LocAxislabel{i})
    xlim(Locylim(i,:))
    ylabel('Zero-crossing')
    set(gca,'Tickdir','out')
end


%% Plot CS amplitudes:

% dimensions of fluo_all_zscore{i_mouse}{i_TT}: time, trials (all), ROIs
figure('Position', [200 50 500 500],'name',['CS amplitudes - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
mean_CS = nan(num_CSTT,num_mice,num_ROI);

for i_TT = 1:num_CSTT
    for i_mouse = 1:num_mice
        mean_CS(i_TT,i_mouse,:) = squeeze(mean(mean(fluo_all_zscore{i_mouse}{i_TT}(ind_CS,:,:),2)));
    end
    h = subplot(3,ceil(num_CSTT/3),i_TT); hold on 
    plot([0 num_ROI],[0 0],'--','color','k','HandleVisibility','off')
    plot(1:num_ROI,squeeze(mean_CS(i_TT,:,:)),'o','color',[0.6 0.6 0.6],'HandleVisibility','off');
    plot(1:num_ROI,squeeze(mean(mean_CS(i_TT,:,:),2)),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2) ;
    title(['TT' num2str(i_TT)])
    ylabel('z-score')
    ylim([-0.5 1])
    xticks(1:num_ROI);
    xtickangle(45)
    set(gca,'TickDir','out','xticklabels',fiberIDs)
end

%% Get half-width of CS peak and of US peak

CS_meantraces = nan(num_TT,num_mice,numel(ind_CS(1)-20:ind_CS(end)+20*2),num_ROI);
US_meantraces = nan(num_TT,num_mice,numel(ind_US_exp_early(1)-20:ind_US_exp(end)+20*2),num_ROI);

time_CS = time_fluo(ind_CS(1)-20:ind_CS(end)+20*2);
time_US = time_fluo(ind_US_exp_early(1)-20:ind_US_exp(end)+20*2);

N = round(SR_fluo{1});   % 1s
alpha = 5;        % (30-1)/(2*5)
g1 = gausswin(N,alpha)/sum(gausswin(N,alpha));

for i_TT = 1:num_TT
    for i_mouse = 1:num_mice
        if i_TT > 1
        % Select only the trials with reward sizes => 5ul
            trials_pos = [];
            for i_R = 4:size(Rsize_pos{i_mouse},2)
                trials_pos = cat(1,trials_pos,Rsize_pos{i_mouse}{i_TT,i_R});
            end
        else
            trials_pos = Rsize_pos{i_mouse}{i_TT,1};
        end

        CS_meantraces(i_TT,i_mouse,:,:) = squeeze(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_CS(1)-20:ind_CS(end)+20*2,trials_pos,:),2));                     % dim: TT, mouse, time, ROI
         
        if i_TT ~= num_TT
            US_meantraces(i_TT,i_mouse,:,:) = squeeze(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_US_exp_early(1)-20:ind_US_exp(end)+20*2,trials_pos,:),2));       % dim: TT, mouse, time, ROI
        else
            US_meantraces(i_TT,i_mouse,:,:) = squeeze(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_US_unexp(1)-30:ind_US_unexp(end)+49,trials_pos,:),2));   % dim: TT, mouse, time ROI
        end
        for i_ROI = 1:num_ROI
            % Smooth with a Gaussian
            CS_meantraces_sm(i_TT,i_mouse,:,i_ROI) = conv(squeeze(CS_meantraces(i_TT,i_mouse,:,i_ROI)),g1,'same');
            US_meantraces_sm(i_TT,i_mouse,:,i_ROI) = conv(squeeze(US_meantraces(i_TT,i_mouse,:,i_ROI)),g1,'same');
        end
    end
end

% Plotting:
figure('Position', [200 50 1700 500],'name',['CS traces - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
for i_TT = 1:num_CSTT
    for i_ROI = 1:num_ROI
        h = subplot(num_TT,num_ROI,i_ROI+num_ROI*(i_TT-1)); hold on 
        plot(time_CS,squeeze(CS_meantraces(i_TT,:,:,i_ROI)),'color',[0.8 0.8 0.8],'HandleVisibility','off');
        plot(time_CS,squeeze(nanmean(CS_meantraces(i_TT,:,:,i_ROI),2)),'color',[0.8 0.5 0.7],'linewidth',2) ;
        plot(time_CS,squeeze(nanmean(CS_meantraces_sm(i_TT,:,:,i_ROI),2)),'color',[0.6 0.3 0.5],'linewidth',1) ;
        if i_TT == 1
            title(fiberIDs{i_ROI})
        end
        ylabel(['TT' num2str(i_TT) ' (z-score'])
        ylim([-0.5 1])
        if i_TT == num_TT
            xlabel('Time (s)')
        end
        set(gca,'TickDir','out')
    end
end
figure('Position', [200 50 1700 500],'name',['US traces - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
for i_TT = 1:num_TT
    for i_ROI = 1:num_ROI
        h = subplot(num_TT,num_ROI,i_ROI+num_ROI*(i_TT-1)); hold on 
        plot(time_US,squeeze(US_meantraces(i_TT,:,:,i_ROI)),'color',[0.8 0.8 0.8],'HandleVisibility','off');
        plot(time_US,squeeze(nanmean(US_meantraces(i_TT,:,:,i_ROI),2)),'color',[0.8 0.5 0.7],'linewidth',2) ;
        plot(time_US,squeeze(nanmean(US_meantraces_sm(i_TT,:,:,i_ROI),2)),'color',[0.6 0.3 0.5],'linewidth',1) ;
        if i_TT == 1
            title(fiberIDs{i_ROI})
        end
        ylabel(['TT' num2str(i_TT) ' (z-score'])
        ylim([-0.5 1])
        if i_TT == num_TT
            xlabel('Time (s)')
        end
        set(gca,'TickDir','out')
    end
end

%% Quantify the peak and the FWHM

fCS = figure('Position', [200 50 1700 500],'name',['CS traces - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
fUS = figure('Position', [200 50 1700 500],'name',['US traces - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
FWHM_CS = nan(num_TT,num_mice,num_ROI);
FWHM_US = nan(num_TT,num_mice,num_ROI);
FWHM_prom_US = FWHM_US;
FWHM_prom_CS = FWHM_US;
for i_TT = 1:num_TT    
    for i_ROI = 1:num_ROI
        for i_mouse = 1:num_mice
            PeakSig_CS = squeeze(CS_meantraces_sm(i_TT,i_mouse,:,i_ROI));
            [~,~,widths,proms] = findpeaks(PeakSig_CS,time_CS,'SortStr','descend','MinPeakProminence',0.2);
            if ~isempty(widths)
                FWHM_CS(i_TT,i_mouse,i_ROI) = widths(1);
                FWHM_prom_CS(i_TT,i_mouse,i_ROI) = proms(1);
            end
            
            PeakSig_US = squeeze(US_meantraces_sm(i_TT,i_mouse,:,i_ROI));
            [~,~,widths,proms] = findpeaks(PeakSig_US,time_US,'SortStr','descend','MinPeakProminence',0.2);
            if ~isempty(widths)
                FWHM_US(i_TT,i_mouse,i_ROI) = widths(1);
                FWHM_prom_US(i_TT,i_mouse,i_ROI) = proms(1);
            end
            
        end
        figure(fCS)
        PeakSig_CS_mean = squeeze(nanmean(CS_meantraces_sm(i_TT,:,:,i_ROI),2));
        subplot(num_TT,num_ROI,i_ROI+num_ROI*(i_TT-1)); hold on
        findpeaks(PeakSig_CS_mean,time_CS,'MinPeakProminence',0.2,'Annotate','extents')
        legend(gca,'off')
        
        figure(fUS)
        PeakSig_US_mean = squeeze(nanmean(US_meantraces_sm(i_TT,:,:,i_ROI),2));
        subplot(num_TT,num_ROI,i_ROI+num_ROI*(i_TT-1)); hold on
        findpeaks(PeakSig_US_mean,time_US,'MinPeakProminence',0.2,'Annotate','extents')
        legend(gca,'off')
    end
end

% Plot FWHM
figure('Position', [200 50 1200 500],'name',['FWHM US - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
for i_TT = 1:num_TT
    h = subplot(3,2,(i_TT-1)*2+1); hold on 
    plot(1:num_ROI,1000*squeeze(FWHM_CS(i_TT,:,:)),'o','color',[0.6 0.6 0.6],'HandleVisibility','off');
    plot(1:num_ROI,1000*squeeze(nanmean(FWHM_CS(i_TT,:,:),2)),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2) ;
    title(['TT' num2str(i_TT) ' CS peak (z-score)'])
    ylabel('FWHM (ms)')
    xticks(1:num_ROI);
    xtickangle(45)
    set(gca,'TickDir','out','xticklabels',fiberIDs)
    
    h = subplot(3,2,(i_TT-1)*2+2); hold on 
    plot(1:num_ROI,1000*squeeze(FWHM_US(i_TT,:,:)),'o','color',[0.6 0.6 0.6],'HandleVisibility','off');
    plot(1:num_ROI,1000*squeeze(nanmean(FWHM_US(i_TT,:,:),2)),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2) ;
    title(['TT' num2str(i_TT) ' US peak (z-score)'])
    ylabel('FWHM (ms)')
    xticks(1:num_ROI);
    xtickangle(45)
    set(gca,'TickDir','out','xticklabels',fiberIDs)
end

% Plot prominences of FWHM
figure('Position', [200 50 1200 500],'name',['FWHM US - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
for i_TT = 1:num_TT
    h = subplot(3,2,(i_TT-1)*2+1); hold on 
    plot(1:num_ROI,squeeze(FWHM_prom_CS(i_TT,:,:)),'o','color',[0.6 0.6 0.6],'HandleVisibility','off');
    plot(1:num_ROI,squeeze(nanmean(FWHM_prom_CS(i_TT,:,:),2)),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2) ;
    title(['TT' num2str(i_TT) ' CS peak (z-score)'])
    ylabel('Prominence (z-score)')
    xticks(1:num_ROI);
    xtickangle(45)
    set(gca,'TickDir','out','xticklabels',fiberIDs)
    
    h = subplot(3,2,(i_TT-1)*2+2); hold on 
    plot(1:num_ROI,squeeze(FWHM_prom_US(i_TT,:,:)),'o','color',[0.6 0.6 0.6],'HandleVisibility','off');
    plot(1:num_ROI,squeeze(nanmean(FWHM_prom_US(i_TT,:,:),2)),'-o','MarkerFaceColor',[1 0.5 0.8],'color',[1 0.5 0.8],'linewidth',2) ;
    title(['TT' num2str(i_TT) ' US peak (z-score)'])
    ylabel('Prominence (z-score)')
    xticks(1:num_ROI);
    xtickangle(45)
    set(gca,'TickDir','out','xticklabels',fiberIDs)
end

% Plot FWHM as a function of location
% CS
figure('Position',[50 50 750 400]);
for i_TT = 1:num_TT
    for i = 1:size(LocAxis,1)
        subplot(num_TT,3,i+(i_TT-1)*3); hold on
        for i_cell = 1:num_ROIs-5   % without TS
            plot(LocAxis(i,i_cell),1000*squeeze(nanmean(FWHM_CS(i_TT,:,i_cell),2)),'o')
        end 
        for i_cell = num_ROIs-4:num_ROIs   % TS only
            plot(LocAxis(i,i_cell),1000*squeeze(nanmean(FWHM_CS(i_TT,:,i_cell),2)),'ok')
        end 
        LocAxis_reg = repmat(LocAxis(i,1:10),num_mice,1);
        FWHM_reg = 1000*squeeze(FWHM_CS(i_TT,:,1:10));
        LocAxis_reg_aux = LocAxis_reg(~isnan(LocAxis_reg));
        FWHM_reg_aux = FWHM_reg(~isnan(LocAxis_reg));
        LocAxis_reg_aux2 = LocAxis_reg_aux(~isnan(FWHM_reg_aux));
        FWHM_reg_aux2 = FWHM_reg_aux(~isnan(FWHM_reg_aux));
        LocAxis_reg = LocAxis_reg_aux2(:);
        FWHM_reg = FWHM_reg_aux2(:);

        [beta, ~, stats] = glmfit(LocAxis_reg,FWHM_reg);
        R = corrcoef(LocAxis_reg,FWHM_reg);
        plot([min(LocAxis_reg)-0.1 max(LocAxis_reg)+0.1], ...
            beta(1) + beta(2)*[min(LocAxis_reg)-0.1 max(LocAxis_reg)+0.1], 'k-')
        title(['R = ' num2str(round(R(1,2),2)) ', p = ' num2str(round(stats.p(2),4))])
        xlabel(LocAxislabel{i})
        xlim(Locylim(i,:))
        ylabel(['TT' num2str(i_TT) ' FWHM CS'])
        set(gca,'Tickdir','out')
    end
end
% US
figure('Position',[50 50 750 400]);
for i_TT = 1:num_TT
    for i = 1:size(LocAxis,1)
        subplot(num_TT,3,i+(i_TT-1)*3); hold on
        for i_cell = 1:num_ROIs-5  % without TS
            plot(LocAxis(i,i_cell),1000*squeeze(nanmean(FWHM_US(i_TT,:,i_cell),2)),'o')
        end 
         for i_cell = num_ROIs-4:num_ROIs  % TS
            plot(LocAxis(i,i_cell),1000*squeeze(nanmean(FWHM_US(i_TT,:,i_cell),2)),'ok')
        end
        LocAxis_reg = repmat(LocAxis(i,1:10),num_mice,1);
        FWHM_reg = 1000*squeeze(FWHM_US(i_TT,:,1:10));
        LocAxis_reg_aux = LocAxis_reg(~isnan(LocAxis_reg));
        FWHM_reg_aux = FWHM_reg(~isnan(LocAxis_reg));
        LocAxis_reg_aux2 = LocAxis_reg_aux(~isnan(FWHM_reg_aux));
        FWHM_reg_aux2 = FWHM_reg_aux(~isnan(FWHM_reg_aux));
        LocAxis_reg = LocAxis_reg_aux2(:);
        FWHM_reg = FWHM_reg_aux2(:);

        [beta, ~, stats] = glmfit(LocAxis_reg,FWHM_reg);
        R = corrcoef(LocAxis_reg,FWHM_reg);
        plot([min(LocAxis_reg)-0.1 max(LocAxis_reg)+0.1], ...
            beta(1) + beta(2)*[min(LocAxis_reg)-0.1 max(LocAxis_reg)+0.1], 'k-')
        title(['R = ' num2str(round(R(1,2),2)) ', p = ' num2str(round(stats.p(2),4))])
        xlabel(LocAxislabel{i})
        xlim(Locylim(i,:))
        ylabel(['TT' num2str(i_TT) ' FWHM US'])
        set(gca,'Tickdir','out')
    end
end

%% This below does not work well

% Normalized each peak to the mean peak in response to uncued water drop of 20 ul
MeanWater20 = nan(num_mice,num_ROI);
for i_mouse = 1:num_mice
    % Mean response to uncued 20ul water delivery:
    R20_pos = Rsize_pos{i_mouse}{num_TT,end};
    MeanWater20(i_mouse,:) = squeeze(nanmean(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_US_unexp,R20_pos,:),2)));
%     MeanWater20(i_mouse,:) = squeeze(nanmean(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_US_unexp,:,:),2)));
end

% Normalize the transients
for i_TT = 1:num_TT
    for i_mouse = 1:num_mice
        CS_meantraces_norm(i_TT,i_mouse,:,:) = squeeze(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_CS(1)-20:ind_CS(end)+20*2,:,:),2))./MeanWater20(i_mouse,:);                     % dim: TT, mouse, time, ROI
        if i_TT ~= num_TT
            US_meantraces_norm(i_TT,i_mouse,:,:) = squeeze(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_US_exp_early(1)-20:ind_US_exp(end)+20*2,:,:),2))./MeanWater20(i_mouse,:);       % dim: TT, mouse, time, ROI
        else
            US_meantraces_norm(i_TT,i_mouse,:,:) = squeeze(nanmean(fluo_all_zscore{i_mouse}{i_TT}(ind_US_unexp(1)-30:ind_US_unexp(end)+49,:,:),2))./MeanWater20(i_mouse,:);   % dim: TT, mouse, time ROI
        end
    end
end
% Plotting:
figure('Position', [200 50 1700 500],'name',['CS traces - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
for i_TT = 1:num_CSTT
    for i_ROI = 1:num_ROI
        h = subplot(num_TT,num_ROI,i_ROI+num_ROI*(i_TT-1)); hold on 
%         plot(time_CS,squeeze(CS_meantraces_norm(i_TT,:,:,i_ROI)),'color',[0.6 0.6 0.6],'HandleVisibility','off');
        plot(time_CS,squeeze(nanmean(CS_meantraces_norm(i_TT,:,:,i_ROI),2)),'color',[0.8 0.5 0.7],'linewidth',2) ;
        if i_TT == 1
            title(fiberIDs{i_ROI})
        end
        ylabel(['TT' num2str(i_TT) ' (z-score'])
%         ylim([-0.5 1])
        if i_TT == num_TT
            xlabel('Time (s)')
        end
        set(gca,'TickDir','out')
    end
end
figure('Position', [200 50 1700 500],'name',['US traces - all TT'],'numbertitle','off', 'Units', 'Pixels');
set(gcf , 'Color', 'w'); hold on
for i_TT = 1:num_TT
    for i_ROI = 1:num_ROI
        h = subplot(num_TT,num_ROI,i_ROI+num_ROI*(i_TT-1)); hold on 
%         plot(time_US,squeeze(US_meantraces_norm(i_TT,:,:,i_ROI)),'color',[0.6 0.6 0.6],'HandleVisibility','off');
        plot(time_US,squeeze(nanmean(US_meantraces_norm(i_TT,:,:,i_ROI),2)),'color',[0.8 0.5 0.7],'linewidth',2) ;
        if i_TT == 1
            title(fiberIDs{i_ROI})
        end
        ylabel(['TT' num2str(i_TT) ' (z-score'])
%         ylim([-0.5 1])
        if i_TT == num_TT
            xlabel('Time (s)')
        end
        set(gca,'TickDir','out')
    end
end


%% Analysis based on concatenation of trials

% Auto-correlation
% See Mohebi et al

x

% Concatenate trials

% To compare CS responses to baseline, concatenate all possible R sizes in the TT
% baselines_compCS = cellfun(@(x) reshape(x,size(x,1),size(x,2),size(x,3)*size(x,4)),baselines,'UniformOutput',false);
% CStrans_compBSL = cellfun(@(x) reshape(x,size(x,1),size(x,2),size(x,3)*size(x,4)),CStrans,'UniformOutput',false);

for i_mouse = 1:num_mice
    
    baselines_compCS = cellfun(@(x) reshape(x,size(x,1),size(x,2),size(x,3)*size(x,4)),baselines,'UniformOutput',false);
    
    for i_TT = 1:num_TT
        
    end
end

% Correlation with movement from the Rotary encoder

end