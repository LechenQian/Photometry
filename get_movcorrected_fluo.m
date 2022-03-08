function fluo_all_movcorr = get_movcorrected_fluo(fluo_all,ind_movcorr,funcChannel,mouse,day,fiberIDs)

num_mice = size(fluo_all,2);
fluo_all_movcorr = cell(1,num_mice);
plotauxfig = 0;

for i_mouse = 1:num_mice
    num_TT = size(fluo_all{i_mouse},2);
    num_ROI = size(fluo_all{i_mouse}{1},3); 
    fluo_all_movcorr{i_mouse} = cell(1,num_TT);
    % For each ROI, get the scatter plot and correlation of all the points
    % in the acquisition:
    h1 = figure('Position', [50 500 2000 300],'name',[mouse{i_mouse} '_' day{i_mouse}],'numbertitle','off', 'Units', 'Pixels');
    h2 = figure('Position', [50 500 2000 300],'name',[mouse{i_mouse} '_' day{i_mouse}],'numbertitle','off', 'Units', 'Pixels');
    
    % Perform the correction at each ROI:
    
% This was just to make a figure for lab meeting:    
%     figure('Position',[675 373 700 601])
%     set(gcf , 'Color', 'w'); hold on
    for i_ROI = 1:num_ROI
        x_isosb_conc = [];
        x_tdT_conc = [];
        y_conc = [];
        for i_TT = 2:num_TT  % Do not use TT1, which only introduces noise, because the mice do not lick
            x_isosb = movmean(fluo_all{i_mouse}{i_TT}(ind_movcorr,:,i_ROI,1),6); 
            x_tdT = movmean(fluo_all{i_mouse}{i_TT}(ind_movcorr,:,i_ROI,3),6); 
            y = movmean(fluo_all{i_mouse}{i_TT}(ind_movcorr,:,i_ROI,funcChannel),6);

            x_isosb = x_isosb(:);
            x_tdT = x_tdT(:);
            y = y(:);
            
            % concatenate data from different trial types
            x_isosb_conc = cat(1,x_isosb_conc,x_isosb);
            x_tdT_conc = cat(1,x_tdT_conc,x_tdT);
            y_conc = cat(1,y_conc,y);
        end
        
        y_conc1 = y_conc;
        y_conc2 = y_conc;
        % Remove points that have amplitude less than 1STD in the
        % functional channel
%         pos_thres_y = y_conc>nanmean(y_conc)+nanstd(y_conc) | y_conc<nanmean(y_conc)-nanstd(y_conc);
%         pos_thres_x1 = x_isosb_conc>nanmean(x_isosb_conc)+nanstd(x_isosb_conc) | x_isosb_conc<nanmean(x_isosb_conc)-nanstd(x_isosb_conc);
%         pos_thres_x2 = x_tdT_conc>nanmean(x_tdT_conc)+nanstd(x_tdT_conc) | x_tdT_conc<nanmean(x_tdT_conc)-nanstd(x_tdT_conc);
%         y_conc1 = y_conc(pos_thres_y & pos_thres_x1);
%         y_conc2 = y_conc(pos_thres_y & pos_thres_x2);
%         x_isosb_conc = x_isosb_conc(pos_thres_y & pos_thres_x1);
%         x_tdT_conc = x_tdT_conc(pos_thres_y & pos_thres_x2);


% This was just to make a figure for lab meeting:
%         plot(x_isosb_conc(~isnan(y_conc))+i_ROI*0.1,'b','linewidth',1)
%         plot(x_tdT_conc(~isnan(y_conc))+i_ROI*0.1,'r','linewidth',1)
%         plot(y_conc(~isnan(y_conc))+i_ROI*0.1,'g','linewidth',1)
%         if i_ROI == num_ROI
%             xlim([17900 18700])
%             set(gca,'Tickdir','out', 'Xcolor', 'w', 'Ycolor', 'w')
%              plot([17910 17910],[1.55 1.55+0.1],'color','k','linewidth',1)
%              text(17912, 1.6,'10% \DeltaF/F_0')
%         end
%     end
        
        % Plotting
        figure(h1)
        subplot(3,num_ROI,i_ROI); hold on
        plot(x_isosb_conc,y_conc1,'.','Markersize',0.2)
        [beta1, ~, stats] = glmfit(x_isosb_conc,y_conc1);
        if stats.p(2) < 0.05
            plot([min(x_isosb_conc)-0.1 max(x_isosb_conc)+0.1], ...
                beta1(1) + beta1(2)*[min(x_isosb_conc)-0.1 max(x_isosb_conc)+0.1], 'k-')
            corr_isosb = 1;
        else 
            corr_isosb = 0;
        end
        R1 = corrcoef([x_isosb_conc,y_conc1],'Rows','complete');
        title(sprintf('%s \n R=%2g',fiberIDs{i_ROI},R1(1,2)));
        if i_ROI == 1
            xlabel('Isosbestic')
            ylabel('Functional')
        end
        
        subplot(3,num_ROI,i_ROI+num_ROI); hold on
        plot(x_tdT_conc,y_conc2,'.','Markersize',0.2)
        [beta2, ~, stats] = glmfit(x_tdT_conc,y_conc2);
        if stats.p(2) < 0.05
            plot([min(x_tdT_conc)-0.1 max(x_tdT_conc)+0.1], ...
                beta2(1) + beta2(2)*[min(x_tdT_conc)-0.1 max(x_tdT_conc)+0.1], 'k-')
            corr_tdT = 1;
        else
            corr_tdT = 0;
        end
        R2 = corrcoef([x_tdT_conc,y_conc2],'Rows','complete');
        title(['R=' num2str(round(R2(1,2),2))])
        if i_ROI == 1
            xlabel('TdTomato')
            ylabel('Functional')
        end
        
        subplot(3,num_ROI,i_ROI+num_ROI*2); hold on
        plot(x_isosb_conc,x_tdT_conc,'.','Markersize',0.2)
        [beta, ~, stats] = glmfit(x_isosb_conc,x_tdT_conc);
        if stats.p(2) < 0.05
            plot([min(x_isosb_conc)-0.1 max(x_isosb_conc)+0.1], ...
                beta(1) + beta(2)*[min(x_tdT_conc)-0.1 max(x_tdT_conc)+0.1], 'k-')
        end
        R = corrcoef([x_isosb_conc,x_tdT_conc],'Rows','complete');
        title(['R=' num2str(round(R(1,2),2))])
        if i_ROI == 1
            xlabel('Isosbestic')
            ylabel('TdTomato')
        end
        
        % Decide whether to use issobestic or tdTomato for correction
        if corr_isosb && ~corr_tdT 
            corrChannel = 1;
            beta_corr = beta1;
        elseif ~corr_isosb && corr_tdT 
            corrChannel = 3;
            beta_corr = beta2;
        elseif corr_isosb && corr_tdT 
            if R1 > R2
               corrChannel = 1;
               beta_corr = beta1;
            else
               corrChannel = 3;
               beta_corr = beta2;
            end
        else
            corrChannel = nan;
            beta_corr = [0 0];
        end
        % Implement correction
        for i_TT = 1:num_TT
            fluo_tocorr = fluo_all{i_mouse}{i_TT}(:,:,i_ROI,funcChannel);
            fluo_ref_isos = fluo_all{i_mouse}{i_TT}(:,:,i_ROI,1);
            fluo_ref_tdT = fluo_all{i_mouse}{i_TT}(:,:,i_ROI,3);
            if ~isnan(corrChannel)
                fluo_ref = fluo_all{i_mouse}{i_TT}(:,:,i_ROI,corrChannel);
                fluo_subtract = (beta_corr(1)+ beta_corr(2)*fluo_ref);
            else
                fluo_ref = 0;
                fluo_subtract = 0;
            end
            fluo_corrected = fluo_tocorr - fluo_subtract;
            fluo_all_movcorr{i_mouse}{i_TT}(:,:,i_ROI) = fluo_corrected;
            
            if i_TT == 2 && plotauxfig 
                figure(h2)
                subplot(3,ceil(num_ROI/3),i_ROI); hold on
                plot(squeeze(mean(fluo_tocorr(:,:,:,:),2)),'k','linewidth',1)
                plot(squeeze(mean(fluo_ref_isos(:,:,:,:),2)),'c')
                plot(squeeze(mean(fluo_ref_tdT(:,:,:,:),2)),'r')
                plot(squeeze(mean(fluo_ref(:,:,:,:),2)),':k')
                plot(squeeze(mean(fluo_subtract(:,:,:,:),2)),'color',[0.5 0 0.8],'linewidth',1)
                plot(squeeze(mean(fluo_corrected(:,:,:,:),2)),'g')
                title(fiberIDs{i_ROI})
                ylabel(['CorrChannel:' num2str(corrChannel)])
                set(gca,'Tickdir','out')
                if i_ROI == num_ROI
                    legend({'GCaMP','Isos','TdT','Ref','ToSub','Corrected'},'NumColumns',2)
                end
            end
        end
    end
end
end