function [Taus, ZCpoints,CellIDTaus,ZCpoints_all] = mat_analysis_OptimismLocation_MSFP(cellid,...
    UStrans,Rsizes,fluotype,protocol,ZCmethod,ASYMmethod,functionType)

[~, ~, ~, CSTT_mostR] = get_protocol_info(protocol);
data = squeeze(UStrans{CSTT_mostR});

%% Colors for asymmetry
ClrmapO = plasma();
ClrmapL = resample(ClrmapO,size(ClrmapO,1)*10,size(ClrmapO,1));
[~,aux2] = max(sum(ClrmapL,2));
aux3 = find(ClrmapL>1,1,'first');
Clrmap = round(ClrmapL(1:min(aux2,aux3),:),4);
Clrmap_asym = linspace(0,1,size(Clrmap,1));

%% Colors for Zero-crossing points
minRsz = min(cellfun(@(x) min(x), Rsizes));
maxRsz = max(cellfun(@(x) max(x), Rsizes));
ClrmapOv = viridis();
ClrmapLv = resample(ClrmapOv,size(ClrmapO,1)*10,size(ClrmapOv,1));
[~,aux2] = max(sum(ClrmapLv,2));
aux3 = find(ClrmapLv>1,1,'first')-1;
Clrmapv = round(ClrmapLv(1:min(aux2,aux3),:),4);
Clrmap_ZC = linspace(minRsz-1,maxRsz+1,size(Clrmapv,1));

figure, hold on
for i = 1:numel(Clrmap_asym)
    plot(Clrmap_asym(i),1,'o','MarkerFaceColor',Clrmap(i,:),'MarkerEdgeColor',Clrmap(i,:))
end
for i = 1:numel(Clrmap_ZC)
    plot(Clrmap_ZC(i),2,'o','MarkerFaceColor',Clrmapv(i,:),'MarkerEdgeColor',Clrmapv(i,:))
end

%% Colors for each session and mouse (color cells based on session and mouse id)
num_cells = numel(cellid);
mice_names = cell(numel(num_cells));
for i_cell = 1:num_cells
    fnToks = strsplit(cellid{i_cell}, '_');
    mice_names{i_cell} = fnToks{1};
    days{i_cell} = fnToks{2};
end
mice_names_uniq = unique(mice_names);
num_mice = length(mice_names_uniq);

possible_CLRmaps = {'Blues','Oranges','Greens','Reds','Greys','Purples'}; 
clrs_cells = nan(num_cells,3);
for i_mouse = 1:num_mice                                                          % one color per mouse
    imouse_sess = unique(days(contains(mice_names,mice_names_uniq{i_mouse})));
    clrs_mouse_sess = brewermap(numel(imouse_sess),possible_CLRmaps{i_mouse});
    for i_mouse_sess = 1:numel(imouse_sess)                                       % one hue per session
        hue_sess = clrs_mouse_sess(i_mouse_sess,:);
        imouse_isess_cells = (contains(cellid,[mice_names_uniq{i_mouse} '_' imouse_sess{i_mouse_sess}]));
        num_cell_clr = sum(imouse_isess_cells);
        clrs_cells(imouse_isess_cells',:) = repmat(hue_sess,[num_cell_clr,1]);                            % all cells of a particular session are the same color and hue
    end
end

Rsz_max = max(cellfun(@(x) max(x),Rsizes));

%% Get Zero crossing points:
[zeroCrossings_Ut, zeroCrossings_Rsz, utilityAxis] = mat_analysis_getZC(data,Rsizes,fluotype,protocol,ZCmethod);
utilityAxisPerCell = utilityAxis - zeroCrossings_Ut;

%% Get asymmetric scaling

[scaleFactNeg, scaleFactNegSE, scaleFactNegN, scaleFactPos,...
    scaleFactPosSE, scaleFactPosN] = get_ScaleFactors(data,...
    utilityAxisPerCell, functionType);

asymM = scaleFactPos ./ (scaleFactNeg + scaleFactPos);

% Limit the values of the scaleFactors, so that they are in the quadrants
% they are supposed to

canEstimate = scaleFactNeg >=0 & scaleFactPos >=0 & ~isnan(scaleFactNeg) & ~isnan(scaleFactPos);

scaleFactNeg_corr = scaleFactNeg;
scaleFactPos_corr = scaleFactPos;
scaleFactNeg_corr(~canEstimate) = nan;
scaleFactPos_corr(~canEstimate) = nan;

asymM_corr = scaleFactPos_corr ./ (scaleFactNeg_corr + scaleFactPos_corr);

% Calculate an asymmetry index based on the angle (not the slope, so that
% it is linear with the difference between the (+) and (-) sides)
% canEstimateT = ~isnan(scaleFactNeg) & ~isnan(scaleFactPos);
% scaleFactNeg_corrT(~canEstimateT) = nan;
% scaleFactPos_corrT(~canEstimateT) = nan;
% scaleFactNeg_corrT = scaleFactNeg; scaleFactNeg_corrT(scaleFactNeg_corrT<0) = 0;
% scaleFactPos_corrT = scaleFactPos; scaleFactPos_corrT(scaleFactPos_corrT<0) = 0;
%  
% asymM_angle_corr = (atan(scaleFactPos_corrT)-atan(scaleFactNeg_corrT))*2/pi;

asymM_angle_corr = (atan(scaleFactPos_corr)-atan(scaleFactNeg_corr))*2/pi;

% Shift this index up so that it is also in [0 1], as the Dabney et al.,
% 2020
asymM_angle_corr = 0.5 + asymM_angle_corr/2;

%% Choose which kind of asym index to use
if strcmp(ASYMmethod,'theta')
    asymM_corr = asymM_angle_corr;
    asymlabel = '$\tau_{\theta} = \frac{(\theta^+-\theta^-)*2}{\pi}\,(*)$';
else
    asymlabel = '$\tau = \frac{\alpha^+}{\alpha^+ + \alpha^-}$';
end

%% Sort according to asymmetric scaling
[asymM_corr_sort,pos_srtAsym] = sort(asymM_corr);
clrs_cells_sortAsym = clrs_cells(pos_srtAsym,:);
zeroCrossings_Rsz_sortAsym = zeroCrossings_Rsz(pos_srtAsym);
cellid_fit = cellid(pos_srtAsym);
num_cells_fit = find(~isnan(asymM_corr_sort),1,'last');

cellid_fit = cellid_fit(1:num_cells_fit);
asymM_corr_sort = asymM_corr_sort(1:num_cells_fit);
zeroCrossings_Rsz_sortAsym = zeroCrossings_Rsz_sortAsym(1:num_cells_fit);
clrs_cells_sortAsym = clrs_cells_sortAsym(1:num_cells_fit,:);

%% Save values to run decoding analysis on individual neurons

Taus = asymM_corr_sort;
ZCpoints = zeroCrossings_Rsz_sortAsym;
CellIDTaus = cellid_fit;
ZCpoints_all = zeroCrossings_Rsz;

%% Plot the correlation of the two

figure; hold on
for i_cell = 1:num_cells_fit
    scatter(asymM_corr_sort(i_cell),zeroCrossings_Rsz_sortAsym(i_cell), 60, 'MarkerFaceColor', clrs_cells_sortAsym(i_cell,:), 'MarkerEdgeColor', 'k')
end
% Fit linear function
[beta, ~, stats] = glmfit(asymM_corr_sort, zeroCrossings_Rsz_sortAsym);
SSRlinear = sum(stats.resid.^2);
R = corr(asymM_corr_sort, zeroCrossings_Rsz_sortAsym, 'rows', 'complete');
plot([-0.05 1.05], beta(1) + beta(2)*[-0.05 1.05], 'k--')
text(-0.05,9,['Linear: RSS=' num2str(round(SSRlinear,1)) ',R=' num2str(round(R,2)) ',p=' num2str(stats.p(2))])
text(-0.05,8.5,['y=' num2str(beta(1)) '+' num2str(beta(2)) 'x'])
% Fit logistic function
fun = @(b,x)b(1)+b(2)*log(x./(1-x));
[coeff,resnorm, residual, ~,~] = lsqcurvefit(fun,[0.1 min(zeroCrossings_Rsz_sortAsym)+1],asymM_corr_sort, zeroCrossings_Rsz_sortAsym);
SSR = sum(residual.^2);
X = 0.0001:0.001:1;
plot(X,coeff(1)+coeff(2)*log(X./(1-X)),'k-')
text(-0.05,8,['Logit: RSS = ' num2str(round(SSR,1))])
text(-0.05,7.5,['y=' num2str(coeff(1)) '+' num2str(coeff(2)) '*log(x/(1-x))'])

% Continue adding info to the plot
text(0.3,Rsz_max+.5,['N = ' num2str(num_cells_fit) ' fit /' num2str(length(asymM)) ' active cells'])
xlim([-0.1 1.1])
ylim([0 Rsz_max+1])
xlabel(asymlabel,'Interpreter','latex')
ylabel('Zero-crossing')
set(gca,'Tickdir','out')
% Make legend for mouse id
text(0.3,Rsz_max,fluotype)
ylim([0 12])


%% Get location of the cells for which asymetry was calculated

% Coordinates of the bottom of the leftmost shank (channels facing us)
MCoor = table_MouseImplantFibers_MSFP;
for i_cell = 1:num_cells_fit
    roinameparts = strsplit(cellid_fit{i_cell},'_');
    ROI = roinameparts{end}; 
    AP(i_cell) = table2array(MCoor(ROI,'AP_target'));              % (1:end-4) is just to remove the ' .mat '
    ML(i_cell) = table2array(MCoor(ROI,'ML_target'));
    DV(i_cell) = table2array(MCoor(ROI,'DV_target'));
%     hemisph = char(table2array(ROI,'ImplantHemisphere'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Similar to the 2P

figure('Position',[50 50 750 400]);
LocAxis = cat(1,AP,ML,DV);
LocAxislabel = {'AP (mm)','ML (mm)', 'DV (mm)'};
Locylim = [-2 2; 0 4; -5 -1];

for i = 1:size(LocAxis,1)
    subplot(2,3,i); hold on
    for i_cell = 1:num_cells_fit
        scatter(LocAxis(i,i_cell),asymM_corr_sort(i_cell), 60, 'MarkerFaceColor', clrs_cells_sortAsym(i_cell,:), 'MarkerEdgeColor', 'k')
        text(LocAxis(i,i_cell),asymM_corr_sort(i_cell),cellid_fit{i_cell}(end-3:end),'Fontsize',5)
    end
    [beta, ~, stats] = glmfit(LocAxis(i,:),asymM_corr_sort);
    R = corr(LocAxis(i,:)',asymM_corr_sort);
    plot([min(LocAxis(i,:))-0.1 max(LocAxis(i,:))+0.1], ...
        beta(1) + beta(2)*[min(LocAxis(i,:))-0.1 max(LocAxis(i,:))+0.1], 'k-')
    title(['R = ' num2str(round(R,2)) ', p = ' num2str(round(stats.p(2),4))])
    xlabel(LocAxislabel{i})
    xlim(Locylim(i,:))
    ylim([-0.1 1.1])
    ylabel(asymlabel,'Interpreter','latex')
end

%% Plot ZC as a function of each axis

for i = 1:size(LocAxis,1)
    subplot(2,3,3+i); hold on
    for i_cell = 1:num_cells_fit
        scatter(LocAxis(i,i_cell),zeroCrossings_Rsz_sortAsym(i_cell), 60, 'MarkerFaceColor', clrs_cells_sortAsym(i_cell,:), 'MarkerEdgeColor', 'k')
        text(LocAxis(i,i_cell),zeroCrossings_Rsz_sortAsym(i_cell),cellid_fit{i_cell}(end-3:end),'Fontsize',5)
    end
    [beta, ~, stats] = glmfit(LocAxis(i,:),zeroCrossings_Rsz_sortAsym);
    R = corr(LocAxis(i,:)',zeroCrossings_Rsz_sortAsym);
    plot([min(LocAxis(i,:))-0.1 max(LocAxis(i,:))+0.1], ...
        beta(1) + beta(2)*[min(LocAxis(i,:))-0.1 max(LocAxis(i,:))+0.1], 'k-')
    title(['R = ' num2str(round(R,2)) ', p = ' num2str(round(stats.p(2),4))])
    xlabel(LocAxislabel{i})
    xlim(Locylim(i,:))
    ylim([minRsz maxRsz])
    ylabel('Zero-crossing')
end

%% Set color according to asymmetric scaling for 3D plot

figure,
subplot(2,2,1); hold on
for i_cell = 1:num_cells_fit
    plot3(AP(i_cell),ML(i_cell),DV(i_cell),'o','MarkerFaceColor',...
        clrs_cells_sortAsym(i_cell,:),'MarkerEdgeColor','k')
    text(AP(i_cell),ML(i_cell),DV(i_cell),cellid_fit{i_cell}(end-3:end),'Fontsize',5)
end
title('Location. Color: mouseID')
xlabel('AP (mm)')
ylabel('ML (mm)')
zlabel('DV (mm)')
grid on
view(74.1032, 19.5947)
title('Color: mouseID')


subplot(2,2,2); hold on
for i_cell = 1:num_cells_fit
    % Set color for each cell based on asymmetry
    [~,ind] = min(abs(Clrmap_asym-asymM_corr_sort(i_cell)));
    asymColor(i_cell,:) = Clrmap(ind,:);
    
    plot3(AP(i_cell),ML(i_cell),DV(i_cell),'o','MarkerFaceColor',...
        asymColor(i_cell,:),'MarkerEdgeColor','k')
    text(AP(i_cell),ML(i_cell),DV(i_cell),['\tau=' num2str(round(asymM_corr_sort(i_cell),1))],'Fontsize',5)
end
title('Location. Color: \tau')
xlabel('AP (mm)')
ylabel('ML (mm)')
zlabel('DV (mm)')
grid on
view(74.1032, 19.5947)
title('Color based on Asym')


%% Set color according to zero-crossing point (all cells, that are not excluded from calculating asymmetric scaling)

[ZC_sort,pos_srt] = sort(zeroCrossings_Rsz);
clrs_cells_sortZC = clrs_cells(pos_srt,:);
cellid_sortZC = cellid(pos_srt);

% Get location for all cells
DV = nan(1,num_cells);
AP = nan(1,num_cells);
ML = nan(1,num_cells);
MCoor = table_MouseImplantFibers_MSFP;
for i_cell = 1:num_cells
    roinameparts = strsplit(cellid{i_cell},'_');
    ROI = roinameparts{end}; 
    AP(i_cell) = table2array(MCoor(ROI,'AP_target'));              % (1:end-4) is just to remove the ' .mat '
    ML(i_cell) = table2array(MCoor(ROI,'ML_target'));
    DV(i_cell) = table2array(MCoor(ROI,'DV_target'));
%     hemisph = char(table2array(ROI,'ImplantHemisphere'));
end

subplot(2,2,3); hold on
for i_cell = 1:num_cells
    plot3(AP(i_cell),ML(i_cell),DV(i_cell),'o','MarkerFaceColor',...
        clrs_cells_sortZC(i_cell,:),'MarkerEdgeColor','k')
    text(AP(i_cell),ML(i_cell),DV(i_cell),cellid_sortZC{i_cell}(end-3:end),'Fontsize',5)
end
title('Location. Color: mouseID')
xlabel('AP (mm)')
ylabel('ML (mm)')
zlabel('DV (mm)')
grid on
view(74.1032, 19.5947)
title('Color: mouseID')

subplot(2,2,4); hold on
for i_cell = 1:num_cells
    % Set color for each cell based on ZC
    [~,ind] = min(abs(Clrmap_ZC-ZC_sort(i_cell)));
    ZCColor(i_cell,:) = Clrmapv(ind,:);
    
    plot3(AP(i_cell),ML(i_cell),DV(i_cell),'o','MarkerFaceColor',...
        ZCColor(i_cell,:),'MarkerEdgeColor','k')
    text(AP(i_cell),ML(i_cell),DV(i_cell),['zc' num2str(round(ZC_sort(i_cell),1))],'Fontsize',5)
end
title('Location. Color: \tau')
xlabel('AP (mm)')
ylabel('ML (mm)')
zlabel('DV (mm)')
grid on
view(74.1032, 19.5947)
title('Color: ZC')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get location of all the cells
% Get location for all cells
Neu_DV = nan(1,num_cells);
Neu_AP = nan(1,num_cells);
Neu_ML = nan(1,num_cells);
MCoor = table_MouseImplantFibers_MSFP;
for i_cell = 1:num_cells
    roinameparts = strsplit(cellid{i_cell},'_');
    ROI = roinameparts{end}; 
    Neu_AP(i_cell) = table2array(MCoor(ROI,'AP_target'));              % (1:end-4) is just to remove the ' .mat '
    Neu_ML(i_cell) = table2array(MCoor(ROI,'ML_target'));
    Neu_DV(i_cell) = table2array(MCoor(ROI,'DV_target'));
%     hemisph = char(table2array(ROI,'ImplantHemisphere'));
end

%% Sort cell coordinates according to ZC
[ZC_sort,pos_srtZC] = sort(zeroCrossings_Rsz);
cellid_sortZC = cellid(pos_srtZC);

Neu_AP_sortZC = Neu_AP(pos_srtZC);
Neu_ML_sortZC = Neu_ML(pos_srtZC);
Neu_DV_sortZC = Neu_DV(pos_srtZC);

% Define colors for sorted ZC
ZCColor = nan(num_cells,3);
for i_cell = 1:num_cells
    [~,ind] = min(abs(Clrmap_ZC-ZC_sort(i_cell)));
    ZCColor(i_cell,:) = Clrmapv(ind,:);    
end

%% Sort cell coordinates according to Asym
cellid_sortAsym = cellid(pos_srtAsym);

Neu_AP_sortAsym = Neu_AP(pos_srtAsym);
Neu_AP_sortAsym = Neu_AP_sortAsym(1:num_cells_fit);
Neu_ML_sortAsym = Neu_ML(pos_srtAsym);
Neu_ML_sortAsym = Neu_ML_sortAsym(1:num_cells_fit);
Neu_DV_sortAsym = Neu_DV(pos_srtAsym);
Neu_DV_sortAsym = Neu_DV_sortAsym(1:num_cells_fit);

% Define colors for sorted Asym
AsymColor = nan(num_cells_fit,3);
for i_cell = 1:num_cells_fit
    [~,ind] = min(abs(Clrmap_asym-asymM_corr_sort(i_cell)));
    AsymColor(i_cell,:) = Clrmap(ind,:);    
end


%% Arranje coordinates for the next plots

LocAxis = cat(1,Neu_AP,Neu_ML,Neu_DV);
LocAxis_sortZC = cat(1,Neu_AP_sortZC,Neu_ML_sortZC,Neu_DV_sortZC);
LocAxis_sortAsym = cat(1,Neu_AP_sortAsym,Neu_ML_sortAsym,Neu_DV_sortAsym);
LocAxislabel = {'AP (mm)','ML (mm)', 'DV (mm)'};
Locylim = [-4.4 -2.6; 0 1.5; -5.2 -3.8];

%% Plot ZC and Asym as a function of each axis

% Asymmetry
figure('Position',[50 50 750 400]);
for i = 1:size(LocAxis_sortAsym,1)
    subplot(2,3,i); hold on
    for i_cell = 1:num_cells_fit
        scatter(LocAxis_sortAsym(i,i_cell),asymM_corr_sort(i_cell), 60, 'MarkerFaceColor', AsymColor(i_cell,:), 'MarkerEdgeColor', 'k')
    end
    [beta, ~, stats] = glmfit(LocAxis_sortAsym(i,:),asymM_corr_sort);
    R = corr(LocAxis_sortAsym(i,:)',asymM_corr_sort);
    plot([min(LocAxis_sortAsym(i,:))-0.1 max(LocAxis_sortAsym(i,:))+0.1], ...
        beta(1) + beta(2)*[min(LocAxis_sortAsym(i,:))-0.1 max(LocAxis_sortAsym(i,:))+0.1], 'k-')
    title(['R = ' num2str(round(R,2)) ', p = ' num2str(round(stats.p(2),4))])
    xlabel(LocAxislabel{i})
    xlim(Locylim(i,:))
    ylim([-0.1 1.1])
    ylabel(asymlabel,'Interpreter','latex')
end

% ZC
for i = 1:size(LocAxis_sortZC,1)
    subplot(2,3,3+i); hold on
    for i_cell = 1:num_cells
        scatter(LocAxis_sortZC(i,i_cell),ZC_sort(i_cell), 60, 'MarkerFaceColor', ZCColor(i_cell,:), 'MarkerEdgeColor', 'k')
    end
    [beta, ~, stats] = glmfit(LocAxis_sortZC(i,:),ZC_sort);
    R = corr(LocAxis_sortZC(i,:)',ZC_sort);
    plot([min(LocAxis_sortZC(i,:))-0.1 max(LocAxis_sortZC(i,:))+0.1], ...
        beta(1) + beta(2)*[min(LocAxis_sortZC(i,:))-0.1 max(LocAxis_sortZC(i,:))+0.1], 'k-')
    title(['R = ' num2str(round(R,2)) ', p = ' num2str(round(stats.p(2),4))])
    xlabel(LocAxislabel{i})
    xlim(Locylim(i,:))
    ylim([minRsz maxRsz])
    ylabel('Zero-crossing (\mul)','Interpreter','latex')
end

%% 3D plots

% Session and mouse color
figure,
subplot(2,2,1); hold on
for i_cell = 1:num_cells
    plot3(Neu_AP(i_cell),Neu_ML(i_cell),Neu_DV(i_cell),'o','MarkerFaceColor',...
        clrs_cells(i_cell,:),'MarkerEdgeColor','k')
    text(Neu_AP(i_cell),Neu_ML(i_cell),Neu_DV(i_cell),cellid{i_cell}(end-3:end),'Fontsize',5)
end
title('Location. Color: mouseID')
xlabel('AP (mm)')
ylabel('ML (mm)')
zlabel('DV (mm)')
grid on
view(74.1032, 19.5947)
title('Location. Color: SessionID')

% Zero-crossing
subplot(2,2,2); hold on
for i_cell = 1:num_cells    
    plot3(Neu_AP_sortZC(i_cell),Neu_ML_sortZC(i_cell),Neu_DV_sortZC(i_cell),'o','MarkerFaceColor',...
        ZCColor(i_cell,:),'MarkerEdgeColor',ZCColor(i_cell,:))
    text(Neu_AP_sortZC(i_cell),Neu_ML_sortZC(i_cell),Neu_DV_sortZC(i_cell),cellid_sortZC{i_cell}(end-3:end),'Fontsize',5)
end
xlabel('AP (mm)')
ylabel('ML (mm)')
zlabel('DV (mm)')
grid on
view(74.1032, 19.5947)
title('Location. Color: ZC')

% Asymmetry
subplot(2,2,3); hold on
for i_cell = 1:num_cells_fit
    plot3(Neu_AP_sortAsym(i_cell),Neu_ML_sortAsym(i_cell),...
        Neu_DV_sortAsym(i_cell),'o','MarkerFaceColor',...
        AsymColor(i_cell,:),'MarkerEdgeColor',AsymColor(i_cell,:))
    text(Neu_AP_sortAsym(i_cell),Neu_ML_sortAsym(i_cell),Neu_DV_sortAsym(i_cell),cellid_sortAsym{i_cell}(end-3:end),'Fontsize',5)
end
xlabel('AP (mm)')
ylabel('ML (mm)')
zlabel('DV (mm)')
grid on
view(74.1032, 19.5947)
title('Location. Color: \tau')
end