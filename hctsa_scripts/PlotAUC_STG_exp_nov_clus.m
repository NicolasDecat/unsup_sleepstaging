%%
clear all
close all

% load('/Users/tand0009/Data/HCTSA_LD/FigJasmine/RawData.mat')
% load('/Users/tand0009/Data/HCTSA_LD/FigJasmine/Psychophysics accuracy CCSHS.mat')
% load('/Users/tand0009/Data/HCTSA_LD/PSG/behav/type1auc_algorithm_table.mat')
load('/Users/nico/Documents/HCTSA/Analysis/AUC/type1auc_human_scorers_table.mat')
load('/Users/nico/Documents/HCTSA/Analysis/AUC/Table_AUC_10iter_v3.mat')

path_raincloud='/Users/nico/Documents/MATLAB/hctsa-master/RainCloudPlots-master/';
addpath(genpath(path_raincloud));
path_export='/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master';
addpath(genpath(path_export));

path_LSCPtools='/Users/nico/Documents/MATLAB/hctsa-master/LSCPtools-master/';
addpath(genpath(path_LSCPtools));

[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);
cl(3, :) = cb(2, :);

%% 

% To adapt my v2 data
SummaryTable.Dataset = num2str(SummaryTable.Dataset);
for b = 1:150
    SummaryTable.Dataset(b,:) = '001';
end
for bb = 151:300
    SummaryTable.Dataset(bb,:) = '005';
end
   
type1auc_algorithm_table=SummaryTable;

Stage = {'Wake' 'N1' 'N2' 'N3' 'REM'};
STG_algo = [0 1 2 3 5];
STG_hum = [1 2 3 4 5];

for S = 1:length(Stage)

stg = Stage{S};

% read into cell array of the appropriate dimensions
data=[];
cond=[1 3];
for i = 1:2
    
    data{i, 3} = [];
    for nset=[1 5]     % D1 ch1, D2 ch1, D1 ch3, D2 ch3
        temp=double(type1auc_algorithm_table.AUC(type1auc_algorithm_table.NumChannels == cond(i) & ...   % All AUC for channel 1 (then channel 3)
            ismember(str2num(type1auc_algorithm_table.Dataset),nset) & ...   % All AUC for Dataset 1 (then D5)   --> AUC of Dataset 001, 1 channel (50 rows, 10 iterations)
            type1auc_algorithm_table.Sleep_stage == STG_algo(S)));                     % Stage Wake (0) only
        tempS=(type1auc_algorithm_table.Dataset((double(type1auc_algorithm_table.NumChannels) == cond(i) & ...  % All dataset for channel 1 (then channel 3)
            ismember(str2num(type1auc_algorithm_table.Dataset),nset) & ...
            type1auc_algorithm_table.Sleep_stage == STG_algo(S)),1:3));
        tempS = cellstr(tempS);
        uS=unique(tempS);
        temp2=[];
        for nSt=1:length(unique(uS))
            temp2=[temp2 nanmean(temp(strcmp(tempS,uS(nSt))))];   % mean AUC (from temp), of all AUC which are from dataset uS (of all temp, in the end)
        end
        data{i, 3} = [data{i, 3} temp2];    %     first row, 3rd col = [mean AUC D1ch1, mean AUC D5ch1].     second row, 3rd col = [mean AUC D1ch2, mean AUC D5ch2]
    end
    
     temp=double(type1auc_algorithm_table.AUC(double(type1auc_algorithm_table.NumChannels) == cond(i) & ...
        ~ismember(str2num(type1auc_algorithm_table.Dataset),[1 5]) & ...   % --> AUC of the 10 other Datasets, 1 channel (50 rows, 10 iterations) 
        type1auc_algorithm_table.Sleep_stage == STG_algo(S)));      
    tempS=(type1auc_algorithm_table.Dataset((double(type1auc_algorithm_table.NumChannels) == cond(i) & ...
        ~ismember(str2num(type1auc_algorithm_table.Dataset),[1 5])& ...
         type1auc_algorithm_table.Sleep_stage == STG_algo(S)),1:3));
    tempS = cellstr(tempS);
    uS=unique(tempS);
    temp2=[];
    for nSt=1:length(unique(uS))
        temp2=[temp2 nanmean(temp(strcmp(tempS,uS(nSt))))];   % mean AUC of all datasets starting with 4, then with 5, 6, 7, 8 (intentional to group them like this? Or shouldn't we make average of each of the 10 datasets
    end
    data{i, 3} = [data{i, 3} temp2];    % at this point: row1, col3, first 2 numbers = AUC for D1 and D5, and the rest 5 is AUC for D4s, D5s, D6s, D7s, D8s (for channel 1). Row2, col3 = same but channel 3
end
for i = 1:2
    data{i, 1} =[];
    dataset{i, 1} =[];
    for nset=[1 5]
        temp=double(type1auc_human_scorers_table.type1auc(double(type1auc_human_scorers_table.channel) == cond(i) & ... 
            type1auc_human_scorers_table.dataset_type=='Novice' & type1auc_human_scorers_table.dataset==num2str(nset) & ...  % All AUC for channel 1 (then channel 3), Novice, and Dataset 1 ((then D5)
            str2double(type1auc_human_scorers_table.sleep_stage) == STG_hum(S)));  % Wake is 1 here
        tempS=(type1auc_human_scorers_table.subject(double(type1auc_human_scorers_table.channel) == cond(i) & ...   
            type1auc_human_scorers_table.dataset_type=='Novice' & type1auc_human_scorers_table.dataset==num2str(nset) & ...   % Get subjects  for AUC rows named above (temp)
            str2double(type1auc_human_scorers_table.sleep_stage) == STG_hum(S)));  
        uS=unique(tempS);
        temp2=[];
        for nSt=1:length(unique(uS))
            temp2=[temp2 nanmean(temp(tempS==uS(nSt)))];   % mean AUC from all subjects involved in dataset 1 (then (), channel 1, Novice (4 subjects)
        end
        data{i, 1} = [data{i, 1} temp2];   % at this point: row1, col3, first 2 numbers = AUC for all novice 4 subjects involved in D1, ch1 (row1, col1: [mean 4 novice AUC D1ch1, mean 4 novice AUC D5ch1]
        dataset{i, 1} = [dataset{i, 1} nset*ones(1,length(temp2))];   
    end
    
    data{i, 2} = [];     % Same but for experts here (only 2, compared to 4 novices)
      dataset{i, 2} =[];
  for nset=[1 5]
    temp=double(type1auc_human_scorers_table.type1auc(double(type1auc_human_scorers_table.channel) == cond(i) & ...
        type1auc_human_scorers_table.dataset_type=='Expert' & type1auc_human_scorers_table.dataset==num2str(nset) & ...  % All AUC for channel 1 (then channel 3), Novice, and Dataset 1 ((then D5)
        str2double(type1auc_human_scorers_table.sleep_stage) == STG_hum(S)));  % Wake is 1 here
    tempS=(type1auc_human_scorers_table.subject(double(type1auc_human_scorers_table.channel) == cond(i) & ...
        type1auc_human_scorers_table.dataset_type=='Expert' & type1auc_human_scorers_table.dataset==num2str(nset)  & ...   % Get subjects  for AUC rows named above (temp)
            str2double(type1auc_human_scorers_table.sleep_stage) == STG_hum(S)));  
    uS=unique(tempS);
    temp2=[];
    for nSt=1:length(unique(uS))
        temp2=[temp2 nanmean(temp(tempS==uS(nSt)))];
    end
         data{i, 2} = [data{i, 2} temp2];    % at this point: row1, col3, first 2 numbers = AUC for all novice 4 subjects involved in D1, ch1 (row1, col1: [mean 4 experts AUC D1ch1, mean 2 experts AUC D5ch1].   For ch3: [mean 3 experts AUC D1ch3, mean 2 experts AUC D5ch3]
         dataset{i, 2} = [dataset{i, 2} nset*ones(1,length(temp2))];
 end
end


figure; format_fig;
set(gcf,'Position',[1   117   400   680]);
for ncond=1:2
    subplot(2, 1, ncond); format_fig;
    
%     h1 = raincloud_plot(data{ncond,1}, 'box_on', 1, 'color', cb(7,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0,'band_width',.04);
%     scatter(h1{2}.XData(dataset{ncond,1}==1),h1{2}.YData(dataset{ncond,1}==1),'Marker','o','MarkerFaceColor',cb(7,:),'MarkerEdgeColor',[1 1 1]*0.5,'LineWidth',1,'SizeData',100,'MarkerFaceAlpha',0.7);
%     scatter(h1{2}.XData(dataset{ncond,1}==5),h1{2}.YData(dataset{ncond,1}==5),'Marker','d','MarkerFaceColor',cb(7,:),'MarkerEdgeColor',[1 1 1]*0.5,'LineWidth',1,'SizeData',100,'MarkerFaceAlpha',0.7);
%     
%     h2 = raincloud_plot(data{ncond,2}, 'box_on', 1, 'color', cb(6,:), 'alpha', 0.5,...
%         'box_dodge', 1, 'box_dodge_amount', .35, 'dot_dodge_amount', .35, 'box_col_match', 0,'band_width',.04);
%     scatter(h2{2}.XData(dataset{ncond,2}==1),h2{2}.YData(dataset{ncond,2}==1),'Marker','o','MarkerFaceColor',cb(6,:),'MarkerEdgeColor',[1 1 1]*0.5,'LineWidth',1,'SizeData',100,'MarkerFaceAlpha',0.7);
%     scatter(h2{2}.XData(dataset{ncond,2}==5),h2{2}.YData(dataset{ncond,2}==5),'Marker','d','MarkerFaceColor',cb(6,:),'MarkerEdgeColor',[1 1 1]*0.5,'LineWidth',1,'SizeData',100,'MarkerFaceAlpha',0.7);
    
    h3 = raincloud_plot(data{ncond,3}, 'box_on', 1, 'color', cb(5,:), 'alpha', 0.5,...
        'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55,...
        'box_col_match', 0,'band_width',.04);
    scatter(h3{2}.XData(3:end),h3{2}.YData(3:end),'SizeData',72,'Marker','s','MarkerFaceColor',cb(5,:),'MarkerEdgeColor',[1 1 1]*0.5,'LineWidth',2,'SizeData',100,'MarkerFaceAlpha',0.7);
%     scatter(h3{2}.XData(1),h3{2}.YData(1),'SizeData',72,'Marker','o','MarkerFaceColor',cb(5,:),'MarkerEdgeColor','r','LineWidth',2,'SizeData',100,'MarkerFaceAlpha',0.7);
%     scatter(h3{2}.XData(2),h3{2}.YData(2),'SizeData',72,'Marker','d','MarkerFaceColor',cb(5,:),'MarkerEdgeColor','r','LineWidth',2,'SizeData',100,'MarkerFaceAlpha',0.7);

    set(gca,'XLim', [0.35 1], 'YLim', ylim.*[.8 1.6]);
    line([1 1]*.5,ylim,'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2)

%     set(h1{2},'SizeData',54)
%     set(h2{2},'SizeData',54)
    set(h3{2},'SizeData',54)
    % xlabel('Type 1 - AUC')
    box off
    format_fig;
    
    % Non visible graph to insert additional legend
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'ok','MarkerEdgeColor','k','LineWidth',2);
    h(2) = plot(NaN,NaN,'dk','MarkerEdgeColor','k','LineWidth',2);
    h(3) = plot(NaN,NaN,'sk','MarkerEdgeColor','k','LineWidth',2);
    
    if ncond == 1
    legend([h3{1}], {'Cluster'},'Location','northwest','FontSize',10);
    title(sprintf('%s',stg))
    end
  
    
    % xlabel('Type 1 - AUC')
    if ncond == 1
        % ylabel({'EEG'})
    else 
        % ylabel({'EEG+EOG+EMG'}) 
        xlabel('Type 1 - AUC')
    end
    
   
end

fpath = '/Users/nico/Documents/HCTSA/Analysis/AUC/AUC_figures/v3/cluster_alone';
saveas(gca,fullfile(fpath,sprintf('CEN_%s stage',stg)),'jpeg')

end

% saveas(gcf, 'thisfigure');
% export_fig('thisfigure.png','-r 30')

% x = export_fig(['/Users/nico/Documents/HCTSA/Analysis/AUC'],'-r 300')
% export_fig(['/Users/tand0009/Work/Documents/Articles/InPrep/HCTSA_LD/Figures' filesep 'Fig_CCSHSclustering_accuracy.fig'])

