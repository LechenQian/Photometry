function plot_pcolorplot_percell(taskPhase,mouse,day,trials,time,means,stds,colors,conditions,CSin,CSout,TraceEnd,revtrial,SR,saraclrmap,num_CSTT,cellids)
dbstop if error
% dimensions of entries in cell 'trials': time,trials,cells
% dimensions of entries in cells 'means' and 'stds': time, cells

filt_time = 0;    % filter across time
filt_trials = 0;  % filter across trials

num_cells = size(trials{1},3);
num_TT = length(trials);
if TraceEnd < 0 
    TraceEnd = CSin+3;
end

maximo = cellfun(@(x) squeeze(nanmax(nanmax(x))),trials,'UniformOutput', false);
maximo = cell2mat(maximo);
maxim = nanmax(maximo,[],2);

minimo = cellfun(@(x) squeeze(nanmin(nanmin(x))),trials,'UniformOutput', false);
minimo = cell2mat(minimo);
minim = nanmin(minimo,[],2);

if ~iscell(means)
    means = mat2cell(means);
end

maxmean_cell = cell2mat(cellfun(@(x) squeeze(nanmax(x))',means,'UniformOutput', false));
minmean_cell = cell2mat(cellfun(@(x) squeeze(nanmin(x))',means,'UniformOutput', false));

% This is needed for the 'fill' function 
time_toplot = -1.8:1/SR:8-1/SR;
pos1 = find(time >= time_toplot(1), 1, 'first');
pos2 = find(time <= time_toplot(end), 1,'last');

% To filt fluo across time
N = SR; % 1s
alpha = 4;  % 15/(2*5) = 
g_time = gausswin(N,alpha)/sum(gausswin(N,alpha));  
% To filt across trials
N = 6; % 1s
alpha = 10;  
g_trials = gausswin(N,alpha)/sum(gausswin(N,alpha)); 

for i_cell = 1:num_cells
    
    if ~isnan(minim(i_cell))
        
        if iscell(cellids)
           cell_figname = cellids{i_cell};
        else
           cell_figname = num2str(cellids(i_cell));
        end
    
        hfig = figure('Position', [200 50 1000 500],'name',[mouse '_' day '_c ' cell_figname],'numbertitle','off', 'Units', 'Pixels');
        set(gcf , 'Color', 'w'); hold on
        clrmap = saracolormap(3,minim(i_cell),maxim(i_cell));

        for i_TT = 1:num_TT
            if filt_time && filt_trials
                datatoplot = conv2(g_time,g_trials,squeeze(trials{i_TT}(:,:,i_cell)),'same')'; % this filters along trials too
            elseif filt_time && ~filt_trials
                datatoplot = conv2(squeeze(trials{i_TT}(:,:,i_cell)),g_time,'same')';
            else
                datatoplot = squeeze(trials{i_TT}(:,:,i_cell))';
            end
            
            datatoplot(1,1) = minim(i_cell);
            datatoplot(1,2) = maxim(i_cell);
            h = subplot(3,num_TT,[i_TT i_TT+num_TT]); hold on
            set(h, 'TickDir', 'out','Fontsize',12,'TickLength',2*(get(gca,'TickLength')),'YTick',0:10:size(datatoplot,2))
            title(conditions{i_TT})
            if i_TT ==1
                ylabel('Trial'); 
                xlabel('Time - Odor (s)')
            end
            
            if taskPhase == 0 || taskPhase == 2    % For Ryu's paper
%                 imagesc(time,1:size(trials{i_TT},2),datatoplot,[-max(abs([datatoplot(1,1),datatoplot(1,2)]))+3 max(abs([datatoplot(1,1),datatoplot(1,2)]))-3]);
                imagesc(time,1:size(trials{i_TT},2)-1,datatoplot,[-nanmax(abs([datatoplot(1,1),datatoplot(1,2)]))+3 nanmax(abs([datatoplot(1,1),datatoplot(1,2)]))-3]);
                colormap yellowblue
                set(gca,'YDir','reverse')
                xlim([-1 5])
                if taskPhase == 0
                    yticks(0:10:size(trials{i_TT},2))
                elseif taskPhase == 2
                    yticks(-1:10:size(trials{i_TT},2)-1)
                    plot([-2 8], [revtrial(i_TT) revtrial(i_TT)]-2,'w','linewidth',2.5)
                    if i_TT <= num_CSTT
                        yticklabels(-10:10:size(trials{i_TT},2)-11);
                    end
                end
                if ~isnan(revtrial(i_TT)) || revtrial(i_TT)~=1
                    if taskPhase == 2
                        ylim([revtrial(i_TT)-10.5 revtrial(i_TT)+18])
                    elseif taskPhase == 0
                        ylim([revtrial(i_TT) revtrial(i_TT)+40])
                    end
                end
                plot([CSin CSin],[1 size(trials{i_TT},2)],'w:','linewidth',1.5)
                plot([TraceEnd TraceEnd],[1 size(trials{i_TT},2)],'w:','linewidth',1.5)
                if i_TT == num_TT
                    colorbar
                    xlim([-1 3])
                end
            else
                pcolor(time,1:size(trials{i_TT},2),datatoplot); shading flat
                if saraclrmap; colormap(clrmap);
                else colormap yellowblue
                end
                xlim([-2 6])
                line([CSin CSin],[1 size(trials{i_TT},2)],'color','k')
                line([TraceEnd TraceEnd],[1 size(trials{i_TT},2)],'color','k')
                line([-2 8], [revtrial(i_TT) revtrial(i_TT)],'color','r')
            end
            xlim([-1 5])

            % Plot mean trace of a cell for each trial type
            h = subplot(3,num_TT,i_TT+2*num_TT); hold on
            if isnan(CSout)     % If there is no CS in this trial type
                sh1 = line([CSin CSin],[minim(i_cell) maxim(i_cell)],'color',[0.9 0.9 0.9],'linewidth',2);
            else
                sh1 = fill([CSin  CSout   CSout CSin], [minim(i_cell) minim(i_cell) maxim(i_cell) maxim(i_cell)],[0.9 0.9 0.9],'EdgeColor','none');
            end
            sh2 = line([TraceEnd TraceEnd],[minim(i_cell) maxim(i_cell)],'color',[0.9 0.9 0.9],'linewidth',2);
            
            if isnan(revtrial(i_TT)) %|| revtrial(i_TT)==1
                % fill is not working
%                 fill([time_toplot fliplr(time_toplot)],[(means{i_TT}(pos1:pos2,i_cell)+stds{i_TT}(pos1:pos2,i_cell))' fliplr((means{i_TT}(pos1:pos2,i_cell)-stds{i_TT}(pos1:pos2,i_cell))')],colors(i_TT,:),'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
                plot(time,means{i_TT}(:,i_cell),'color',colors(i_TT,:),'linewidth',1)
                ylim([min(means{i_TT}(:,i_cell))-0.05 max(means{i_TT}(:,i_cell))+0.05])
                set(h, 'TickDir', 'out','Fontsize',12)
                xlabel('Time from CS onset (s)')
                if i_TT > num_CSTT
                    xlim([-1 3])
                else
                    xlim([-1 5])
                end
            else
                % General Figure for analysis
                % Figure for Ryu's papers
                if taskPhase == 2 || taskPhase == 0
                    plot(time(2:end),nanmean(datatoplot(revtrial(i_TT):revtrial(i_TT)+3,(2:end))),'color',[145 34 179]/255,'linewidth',2)
                    plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+4:revtrial(i_TT)+7,(2:end))),'color',[108 70 179]/255,'linewidth',2)
                    plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+8:revtrial(i_TT)+11,(2:end))),'color',[74 105 179]/255,'linewidth',2)
                    plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+12:revtrial(i_TT)+15,(2:end))),'color',[37 142 179]/255,'linewidth',2)
                    plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+16:revtrial(i_TT)+19,(2:end))),'color',[0 179 179]/255,'linewidth',2)
%                 elseif taskPhase == 0
%                     plot(time(2:end),nanmean(datatoplot(revtrial(i_TT):revtrial(i_TT)+8,(2:end))),'color',[145 34 179]/255,'linewidth',2)
%                     plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+9:revtrial(i_TT)+16,(2:end))),'color',[108 70 179]/255,'linewidth',2)
%                     if size(datatoplot,1) >= revtrial(i_TT)+24
%                         plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+17:revtrial(i_TT)+24,(2:end))),'color',[74 105 179]/255,'linewidth',2)
%                     end
%                     if size(datatoplot,1) >= revtrial(i_TT)+32
%                         plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+25:revtrial(i_TT)+32,(2:end))),'color',[37 142 179]/255,'linewidth',2)
%                     end
%                     if size(datatoplot,1) > revtrial(i_TT)+33
%                         plot(time(2:end),nanmean(datatoplot(revtrial(i_TT)+33:end,(2:end))),'color',[0 179 179]/255,'linewidth',2)
%                     end
                else
                    if ~isempty(datatoplot(1:revtrial(i_TT)-1,:))
                        fill([time(pos1+50:pos2-100) fliplr(time(pos1+50:pos2-100))],[nanmean(datatoplot(1:revtrial(i_TT)-1,pos1+50:pos2-100))+nanstd(datatoplot(1:revtrial(i_TT)-1,pos1+50:pos2-100))/length(1:revtrial(i_TT)-1)' fliplr(nanmean(datatoplot(1:revtrial(i_TT)-1,pos1+50:pos2-100))-nanstd(datatoplot(1:revtrial(i_TT)-1,pos1+50:pos2-100))/length(1:revtrial(i_TT)-1)')],max(0,min(1,colors(i_TT,:)+0.15)),'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
                        plot(time,nanmean(datatoplot(1:revtrial(i_TT)-1,:)),'color',max(0,min(1,colors(i_TT,:)+0.15)),'linewidth',1)
                    end
                    erroaux = nanstd(datatoplot(revtrial(i_TT):end,pos1:pos2),[],1)/sqrt(size(datatoplot(revtrial(i_TT):end,pos1:pos2),1));
                    mediaerroaux = nanmean(datatoplot(revtrial(i_TT):end,pos1:pos2));
                    mediaerro = mediaerroaux(~isnan(erroaux));
                    erro = erroaux(~isnan(erroaux));
                    time_toplotaux = time(pos1:pos2);
                    time_toplot = time_toplotaux(~isnan(erroaux));
                    
                    fill([time_toplot fliplr(time_toplot)],[mediaerro+erro fliplr(mediaerro-erro)],colors(i_TT,:),'facealpha',0.2,'EdgeColor','none','FaceVertexAlphaData',0.2)
                    plot(time(2:end),nanmean(datatoplot(revtrial(i_TT):end,(2:end))),'color',max(0,min(1,colors(i_TT,:))),'linewidth',1)
                end  
                if i_TT > num_CSTT
                    xlim([-1 3])
                else
                    xlim([-1 5])
                end
                set(h, 'TickDir','out','TickLength',2*(get(gca,'TickLength')),'Fontsize',12)
%                 ylim([-8 15])
%                 ylim([-1.5 3])
                ylim([nanmin(minmean_cell(i_cell,:))-0.25 nanmax(maxmean_cell(i_cell,:))+0.25])
%                 ylim([-10 20])
                xlabel('Time - Odor (s)')
            end
        end
    end
%     savefolder = 'C:\Users\sara\Dropbox (Uchida Lab)\Manuscripts\Ryu_2021\Figs_20210504\MainFig1_Revers_Reversal_ExampleCells_SM103_20200205/';
%     saveas(gca,[savefolder mouse '_' day '_cell' num2str(i_cell) '_zscore'],'fig')
% saveas(hfig,[hfig.Name '_baselineSub.fig'])
% saveas(hfig,[hfig.Name '_baselineSub.tif'])
end

end