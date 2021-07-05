 

%% Identify the top features

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))
    
    %%% Features reordering;
    means = mean(Per_correct_mean);
    [~,I] = sort((means)','descend');
    Per_correct_mean = Per_correct_mean(:,I);

    %%% Best X features 
    Top_Feat(:,D) = I(1:20);
    
end

gpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Top_50_feat';
save(fullfile(gpath,'Top_Feat.mat'))

%% Features reordering based on TS_CLUSTER (op_clust)


Subs = {'001'};  % '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001'}; % ,'005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))
    
    % Load dataset 1
    load('HCTSA_N.mat','op_clust')
    
    Per_correct_mean = Per_correct_mean(:,(op_clust.ord)');

end


%% Return the top features from averaged matrix (WB feat)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Equivalent indices for WB-features-only data
load('HCTSA.mat', 'Operations')  

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat);   % remove SV features
Operations = Operations(equi_Top_Feat,:);               % get Operations with WB features only

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

% Features reordering: from best to worst feature
means = mean(AveragedMatrix_excl);  
[~,I] = sort((means)','descend');
AveragedMatrix = AveragedMatrix_excl(:,I);

TOP = 40;

Top_Feat = I(1:TOP);   

Top_name = CodeString(Top_Feat,1);
Top_key = Keywords(Top_Feat,1);
Top_ID = YLabel(Top_Feat,1);
Top_mean = mean(AveragedMatrix(:,1:TOP))';

Top_40 = [Top_name Top_key Top_ID num2cell(Top_mean)];


%% Same as above but top features for each classifier

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Get 'Operations of WB features only
load('HCTSA.mat', 'Operations')  

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat);     % remove SV features
Operations = Operations(equi_Top_Feat,:);                 % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(AveragedMatrix_excl(C,:)','descend');
    Top_Feats{C} = I(1:10);   
                     
    % Get the name and keyword associated with these features
    Top_name(1:10) = CodeString(Top_Feats{C},1);
    Top_key(1:10) = Keywords(Top_Feats{C},1);
    Top_ID(1:10) = YLabel(Top_Feats{C},1);
    Top_mean(1:10) = AveragedMatrix_excl(C,Top_Feats{C});

    Top_10{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];

end


%% Get top features EXCLUDING the special values features (all further scripts are now updated, taking into account exclusion of SV features)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

load('HCTSA.mat', 'Operations')   % Simply to get the list of 7749 features 

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat); % get the equivalent ranking of top feat after special-value features removed
Operations = Operations(equi_Top_Feat,:);  % get Operations of well-behaved features only

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

%%% and replace AverageMatrix and Per_correct_mean_D by AverageMatrix_excl and Per_correct_mean_D_excl


%% Line plot of accuracy from all features sorted from best to worst

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
AveragedMatrix = AveragedMatrix_excl;

% Features reordering: from best to worst feature
means = mean(AveragedMatrix);  % AveragedMatrix = 7749-long including special value features
[~,I] = sort((means)','descend');
AveragedMatrix = AveragedMatrix(:,I);

% Get index of 'red' features (features to remove) and remove them
x = find(AveragedMatrix(1,:) == 0);
AveragedMatrix(:,x) = [];

% Get the average over all 10 classifiers, across features
MeanFeat = mean(AveragedMatrix);

% Line plot
figure; ax = gca;
plot(MeanFeat,'LineWidth',2)
xlabel('Features')
ylabel('Singl-feature classification Accuracy (%)')

ax.XTick = 0:500:5603;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:5603,'uni',0);
labels = string(ax.XTickLabels); 
% labels(2:2:end) = NaN;    % remove every other label
ax.XTickLabels = labels; 
ax.FontSize = 12; 
grid on

% % Plot std (variance across classifiers)
% figure; ax = gca;
% stdshade(AveragedMatrix,0.5,'r');
% title('Mean accuracy across features')
% xlabel('features sorted from best to worst')
% ylabel('Percentage Accuracy')
% ax.XTick = 1:50:5603;
% ax.XTickLabels = arrayfun(@(a)num2str(a),0:50:5603,'uni',0);
% labels = string(ax.XTickLabels); 
% labels(2:2:end) = NaN;    % remove every other label
% ax.XTickLabels = labels; 
% grid on
% xlim([0 5650])
% 
% % Plot std (variance across datasets)
% load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
% 
% for D = 1:12
%     Var(D,:) =  Per_correct_mean_D_excl{1,D}(5,:);
%     Var2(D,:) =  Per_correct_mean_D_excl{1,D}(7,:);
%     Var3(D,:) =  Per_correct_mean_D_excl{1,D}(9,:);
% end
% 
% MEAN = mean(Var);
% [~,I] = sort((MEAN)','descend');
% Var = Var(:,I);
% MEAN2 = mean(Var2);
% [~,Y] = sort((MEAN2)','descend');
% Var2 = Var2(:,Y);
% MEAN3 = mean(Var3);
% [~,Z] = sort((MEAN3)','descend');
% Var3 = Var3(:,Z);
% 
% figure; subplot(3,1,1) 
% stdshade(Var,0.5,'b');
% title('N1/N2')
% xlabel('features')
% ylabel('Percentage Accuracy')
% grid on
% xlim([0 5650])
% ylim([0 100])
% 
% subplot(3,1,2) 
% stdshade(Var2,0.5,'g');
% title('N1/REM')
% xlabel('features')
% ylabel('Percentage Accuracy')
% grid on
% xlim([0 5650])
% ylim([0 100])
% 
% subplot(3,1,3) 
% stdshade(Var3,0.5,'r');
% title('N2/REM ')
% xlabel('features')
% ylabel('Percentage Accuracy')
% grid on
% xlim([0 5650])
% ylim([0 100])


% Line plot but zoom in to top 100 only
figure; ax = gca;
plot(MeanFeat(1:120),'LineWidth',2)
title('')
xlabel('')
ylabel('')

ax.XTick = 1:10:110;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:10:110,'uni',0);
labels = string(ax.XTickLabels); 
ax.XTickLabels = labels; 
ax.FontSize = 12; 
grid on

% xlines for top features
xline(11,'--','top 10','LineWidth',2,'Color','#A2142F','FontSize',13);
xline(41,'--','top 40','LineWidth',2,'Color','#A2142F','FontSize',13);
xline(101,'--','top 100','LineWidth',2,'Color','#A2142F','FontSize',13);


% %%%% Percentile g
% means = means(1,I);
% Y = prctile(means,99)  % Take 1% (53 Features) of best features -> Take top 50 ?




%% Plot the top features - Ben's functions

%%%%% bits of code 

% load('HCTSA_N.mat')
% load('ccshs_1800001_annot.mat', 'sleepstage')   
% 
% % Get EEG-only time series (1374 for 001, 1442 for 005)
% TimeSeries = TimeSeries(1:1374,:);
% TS_DataMat = TS_DataMat(1:1374,:);
% TS_Quality = TS_Quality(1:1374,:);
% save HCTSA_N.mat

% To use to convert TimeSeries.Data into cell array (for TS_PlotDataMatrix)
% Data = [];
% for D = 1:size(TimeSeries.Data,1)
%     Data{D,1} = TimeSeries.Data(D,:);
% end
% 
% TimeSeries.Data = Data;

%% on obtaining Per_correct_mean_D / Per_correct_mean_D_SVM / Per_correct_mean_D_excl

% Copy pasted script to get 7749-matrix
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset

InsertCol = zeros(10,1);

Subs = {'001'}; % '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once
% (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_each/Per_correct_mean(Dataset %s)',sub))
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering

    % I'll insert both the commonly and specifically removed features for
    % this dataset
    spec_common = sort([val specifically_removed{1,y(D)}]);
    
%     for x = 1:length(spec_common)   % comment the for loop if want to have Per_correct_mean_D
% 
%         Part1 = Per_correct_mean(:,1:spec_common(x)-1);
%         Part2 = [InsertCol Per_correct_mean(:,spec_common(x):end)];
%         Per_correct_mean = ([Part1 Part2]);
% 
%     end
    
    % If want Per_correct_mean_D_excl
%     load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')  % For each of the 12 cells, only features specifically removed in the corresponding dataset
%     Per_correct_mean(:,spec_and_common_feat) = [];
    Per_correct_mean_D{D} = Per_correct_mean;
    
%   % If want Per_correct_mean_D_excl
%     Per_correct_mean_D_SVM{D} = Per_correct_mean;
    
end

%% SVM 

%%% 1) Section above ("on obtaining Per_correct_mean_D / or
%%% Per_correct_mean_D_SVM") to get Per_correct_mean_D_SVM
%%% 2) Remove SV features:
%%% Per_correct_mean_D_SVM{1,1}(:,spec_and_common_feat) = []; to get
%%% 10x5603 matrix


%% Top 40 features for one Dataset (only WB Features)

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
end

% Per_correct_mean = Per_correct_mean_D_SVM{y};  % if dealing with svm
Per_correct_mean = Per_correct_mean_D_excl{1,y};

% Reorder features: from yielding theg highest to lowest accuracy
means = mean(Per_correct_mean);
[~,I] = sort((means)','descend');
Per_correct_mean = Per_correct_mean(:,I);
    
% Get the best 40 features
Top_Feat = I(1:40);   

% Load Data
load HCTSA_N.mat

% From Operations, take only well behaved features
equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
Operations = Operations(Idx_WB_Feat,:);                      % 'Operations' with WB features only

% Get the name and keyword associated with these features
CodeString = Operations.CodeString;  
Keywords = Operations.Keywords;
YLabel = Operations.ID;

Top_name(1:40) = CodeString(Top_Feat,1);
Top_key(1:40) = Keywords(Top_Feat,1);
Top_ID(1:40) = YLabel(Top_Feat,1);

Top_40 = [Top_name' Top_key' num2cell(Top_ID')];

% Top_mean (mean over classifiers)
for F = 1:length(Top_40)
    for C = 1:10
        Top_40_Acc(F,C) = mean(Per_correct_mean(C,F));
    end
end

Top_mean = mean(Top_40_Acc');
Top_40 = [Top_40 num2cell(Top_mean')];


%%%%%%%% Plot the correlation matrix %%%%%%%%%


% Get pairwise similarity matrix
numTopFeatures = 40;    % Number of top features
op_ind = Top_ID';       % indices of top 40 features (among the 7k), used if HCTSA_N = false

% Compute correlations based on hctsa responses
HCTSA_N = true;   % true if want to plot normalized hctsa values; else plot HCTSA.mat values

if HCTSA_N
    EEGonly = 1:size(TS_DataMat,1)/7;
    TS_DataMat = TS_DataMat(EEGonly,Idx_WB_Feat);      % Only WB features and EEG channels 
    Dij = BF_pdist(TS_DataMat(:,Top_Feat)','abscorr');  
else
    load('HCTSA.mat', 'TS_DataMat')
    EEGonly = 1:size(TS_DataMat,1)/7;
    Dij = BF_pdist(TS_DataMat(EEGonly,op_ind)','abscorr');    
end
    
distanceMetric = 'abscorr';
clusterThreshold = 0.2;      % threshold at which split into clusters

% Ylabels
Top_mean = num2cell(mean(Top_40_Acc'));  
YLabel = [];
for i = 1:length(Top_ID)
    YLabel = [YLabel {sprintf('[%s] %s (%1.1f%%)',Top_key{i},Top_name{i},Top_mean{i})}];
end

% Plot the correlation matrix
[~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',YLabel);
title(sprintf('Dependencies between %u top features (organized into %u clusters)',...
                        numTopFeatures,length(cluster_Groupi)))
                    

%% Corr matrix averaged over all 12 datasets

%%%%% Get the AveragedMatrix with all 12 datasets and their 5603 WB features

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Get 'Operations of WB features only, from HCTSA_N
load HCTSA.mat   
Operations_ID = [Operations.ID].';

equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); % remove SV features
Operations = Operations(equi_Top_Feat,:);                      % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

% Features reordering: from best to worst feature
means = mean(AveragedMatrix_excl);  
[~,I] = sort((means)','descend');
AveragedMatrix = AveragedMatrix_excl(:,I);

% Get best 40 features 
Top_Feat = I(1:40); 

Top_name(1:40) = CodeString(Top_Feat,1);
Top_key(1:40) = Keywords(Top_Feat,1);
Top_ID(1:40) = YLabel(Top_Feat,1);

Top_40 = [Top_name' Top_key' Top_ID'];

% Get mean over classifiers
for F = 1:length(Top_40)
    for C = 1:10
        Top_40_Acc(F,C) = mean(AveragedMatrix(C,F));
    end
end

Top_mean = mean(Top_40_Acc');

Top_40 = [Top_40 num2cell(Top_mean')];


%%%% Now get the hctsa values averaged over all datasets

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Load TS_DataMat 
    load('HCTSA_N.mat','Operations','TS_DataMat')   % Simply to get the list of 7749 features 
    
    % Get equivalent indices (only WB features)
    equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
    Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
    
    % Get DataMat
    EEGonly = 1:size(TS_DataMat,1)/7;
    DataMat{D} = TS_DataMat(EEGonly,Idx_WB_Feat);
    
end


% Trim the TS_DataMat of each file to have their length matched, before
% averagaging (here, shortest file has 1088 time series (EEG only) 
for x = 1:length(DataMat)
    DataMat{1,x} = DataMat{1,x}(1:1088,:);  % Multiply 1088 by 7 to include all 7 channels, and not only EEG
end

% Average the DataMat files
StackedDataMat = cat(3,DataMat{1},DataMat{2},DataMat{3},DataMat{4},DataMat{5},DataMat{6},DataMat{7},DataMat{8},DataMat{9},DataMat{10},DataMat{11},DataMat{12});
AveragedDataMat = mean(StackedDataMat,3);

%%%%%%%% Plot the correlation matrix

% Get pairwise similarity matrix
numTopFeatures = 40;  % Number of top features
op_ind = cell2mat(Top_ID'); % indices of top 40 features (among the 7k)

% Compute correlations based on hctsa responses
HCTSA_N = true;   % true if want to plot normalized hctsa values; else plot HCTSA.mat values

if HCTSA_N
    EEGonly = 1:size(TS_DataMat,1)/7;
    TS_DataMat = TS_DataMat(EEGonly,Idx_WB_Feat);      % Only WB features and EEG channels 
    Dij = BF_pdist(TS_DataMat(:,Top_Feat)','abscorr');  
else
    load('HCTSA.mat', 'TS_DataMat')
    EEGonly = 1:size(TS_DataMat,1)/7;
    Dij = BF_pdist(TS_DataMat(EEGonly,op_ind)','abscorr');    
end

distanceMetric = 'abscorr';
clusterThreshold = 0.2; % threshold at which split into clusters

% Ylabels
Top_mean = num2cell(mean(Top_40_Acc'));  % Change it for later in the Ylabels section
YLabel = [];

for i = 1:length(Top_ID)
    YLabel = [YLabel {sprintf('[%s] %s (%1.1f%%)',Top_key{i},Top_name{i},Top_mean{i})}];
end

% Plot
[~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',YLabel);
title('All classifiers','FontSize',20)

fpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Top_Features/ind_classifier(all_datasets)/figs';
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
%export_fig([fpath filesep 'AllD_allC_unnorm'],'-r 300')

%% Corr matrix on specific classifiers

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'005'};  % '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
end

Per_correct_mean = Per_correct_mean_D_excl{y};

% Load Data
load HCTSA_N.mat

% From Operations, take only well behaved features
equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
Operations = Operations(Idx_WB_Feat,:);  

% Get the name and keyword associated with these features
CodeString = Operations.CodeString;  
Keywords = Operations.Keywords;
YLabel = Operations.ID;

TOP = 1;

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(Per_correct_mean(C,:)','descend');
    Top_Feats{C} = I(1:TOP);   
                     
    % Get the name and keyword associated with these features
    Top_name(1:TOP) = CodeString(Top_Feats{C},1);
    Top_key(1:TOP) = Keywords(Top_Feats{C},1);
    Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
    Top_mean(1:TOP) = Per_correct_mean(C,Top_Feats{C});

    Top_40{C} = [Top_name' Top_key' num2cell(Top_ID') num2cell(Top_mean')];

end

% [7104 6329 142 2065 2673 92 1009 7011];

% Indices of classifiers
wake = 0; N1 = 1; N2 = 3; N3 = 3; rem = 5;
CLASSIFIER = {[wake,N1] [wake,N2] [wake,N3] [wake,rem] [N1,N2] [N1,N3] [N1,rem] [N2,N3] [N2,rem] [N3,rem]};

% Which classifier do we want to plot?
ClNum = 3;  
Classif = Top_40{1,ClNum};  % Get the top features of the chosen classifier
Top_Feat = Top_Feats{1,ClNum};

% hctsa values of this classif
load(sprintf('ccshs_1800%s_annot.mat',sub), 'sleepstage')   
stage1 = find(sleepstage == CLASSIFIER{ClNum}(1));  % indices of epochs that are labeled stage1
stage2 = find(sleepstage == CLASSIFIER{ClNum}(2));  % indices of epochs that are labeled stage2

DataMat = sort([stage1;stage2]);  % indices of epochs stage1 + stage2

%%%%%%%% Plot the correlation matrix

% Get pairwise similarity matrix
numTopFeatures = 40;         % Number of top features
op_ind = Classif(:,3)';     % 3 = top_ID colum: indices of top 40 features (among the 7k)
op_ind = cell2mat(op_ind);

% Compute correlations based on hctsa responses
HCTSA_N = true;   % true if want to plot normalized hctsa values; else plot HCTSA.mat values

if HCTSA_N
    EEGonly = 1:size(TS_DataMat,1)/7;
    TS_DataMat = TS_DataMat(EEGonly,Idx_WB_Feat);      % Only WB features and EEG channels 
    Dij = BF_pdist(TS_DataMat(:,Top_Feat)','abscorr');  
else
    load('HCTSA.mat', 'TS_DataMat')
    EEGonly = 1:size(TS_DataMat,1)/7;
    Dij = BF_pdist(TS_DataMat(EEGonly,op_ind)','abscorr');    
end

distanceMetric = 'abscorr';
clusterThreshold = 0.2; % threshold at which split into clusters

% Ylabels
Top_mean = Classif(:,4)';  
YLabel = [];
for i = 1:length(Top_ID)
    YLabel = [YLabel {sprintf('[%s] %s (%1.1f%%)',Classif{i,2},Classif{i,1},Classif{i,4})}];
end

% Plot the correlation matrix
[~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',YLabel);
title(sprintf('Dependencies between %u top features (organized into %u clusters)',...
                        numTopFeatures,length(cluster_Groupi)))
           
                    
%% Corr matrix averaged over all 12 datasets on specific classifiers

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Get 'Operations of WB features only, from HCTSA_N
load HCTSA.mat   
Operations_ID = [Operations.ID].';

equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); % remove SV features
Operations = Operations(equi_Top_Feat,:);                      % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

TOP = 40;

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(AveragedMatrix_excl(C,:)','descend');
    Top_Feats{C} = I(1:TOP);   
                    
    % Get the name and keyword associated with these features
    Top_name(1:TOP) = CodeString(Top_Feats{C},1);
    Top_key(1:TOP) = Keywords(Top_Feats{C},1);
    Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
    Top_mean(1:TOP) = AveragedMatrix_excl(C,Top_Feats{C});

    Top_40{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];
end

% Indices of classifiers
wake = 0; N1 = 1; N2 = 3; N3 = 3; rem = 5;
CLASSIFIER = {[wake,N1] [wake,N2] [wake,N3] [wake,rem] [N1,N2] [N1,N3] [N1,rem] [N2,N3] [N2,rem] [N3,rem]};


%%%% Now get the hctsa values averaged over all datasets

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Load TS_DataMat 
    load('HCTSA_N.mat','Operations','TS_DataMat')   % Simply to get the list of 7749 features 
    
    % Get equivalent indices (only WB features)
    equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
    Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
    
    % Get DataMat
    EEGonly = 1:size(TS_DataMat,1)/7;
    DataMat{D} = TS_DataMat(EEGonly,Idx_WB_Feat);
    
end

% Trim the TS_DataMat of each file to have their length matched, before
% averagaging (here, shortest file has 1088 time series (EEG only) 
for x = 1:length(DataMat)
    DataMat{1,x} = DataMat{1,x}(1:1088,:);  % Multiply 1088 by 7 to include all 7 channels, and not only EEG
end

% Average the DataMat files
StackedDataMat = cat(3,DataMat{1,1},DataMat{1,2},DataMat{1,3},DataMat{1,4},DataMat{1,5},DataMat{1,6},DataMat{1,7},DataMat{1,8},DataMat{1,9},DataMat{1,10},DataMat{1,11},DataMat{1,12});
AveragedDataMat = mean(StackedDataMat,3);

%%%%%%%% Plot the correlation matrix

for C = 1:10
    
    % Which classifier do we want to plot?
    ClNum = C;  
    Classif = Top_40{1,ClNum};  % Get the top features of the chosen classifier
    Top_Feat = Top_Feats{1,ClNum};

    % Get pairwise similarity matrix
    numTopFeatures = TOP;  % Number of top features
    op_ind = Classif(:,3)'; % 3 = top_ID colum: indices of top 40 features (among the 7k)
    op_ind = cell2mat(op_ind);

    % Compute correlations based on hctsa responses
    TS_DataMat = AveragedDataMat;      % Only WB features and EEG channels 
    Dij = BF_pdist(TS_DataMat(:,Top_Feat)','abscorr');  

    distanceMetric = 'abscorr';
    clusterThreshold = 0.2;      % threshold at which split into clusters

    % Ylabels
    Top_mean = Classif(:,4)';  % Change it for later in the Ylabels section
    YLabel = [];

    for i = 1:length(Top_ID)
        YLabel = [YLabel {sprintf('[%s] %s (%1.1f%%)',Classif{i,2},Classif{i,1},Classif{i,4})}];
    end

    classifier = {'Wake/N1','Wake/N2','Wake/N3','Wake/REM','N1/N2','N1/N3','N1/REM','N2/N3','N2/REM','N3/REM'};  

    % Plot the correlation matrix
    [~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                            'whatDistance',distanceMetric,...
                            'objectLabels',YLabel);
    title(sprintf('%s   Top %s features',classifier{C},string(numTopFeatures)))


    % Save figure
    classifiers = {'Wake_N1','Wake_N2','Wake_N3','Wake_REM','N1_N2','N1_N3','N1_REM','N2_N3','N2_REM','N3_REM'};  
    thisClassif = classifiers(C);


%     % Save
    fpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Top_Features/ind_classifier(all_datasets)/figs';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
  %  export_fig([fpath filesep sprintf('CM_top40_%s(allD)',string(thisClassif))],'-r 300')

    
end

%% Get the top features for each classifier and each dataset
                   
% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
end

Per_correct_mean = Per_correct_mean_D_excl{y};

% Load Data
load HCTSA_N.mat

% From Operations, take only well behaved features
equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
Operations = Operations(Idx_WB_Feat,:);  

% Get the name and keyword associated with these features
CodeString = Operations.CodeString;  
Keywords = Operations.Keywords;
YLabel = Operations.ID;
TOP = 40;

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:12
    
    sub = Subs{D};
   
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))      

    for C = 1:10

        % Features reordering: from best to worst feature
        [~,I] = sort(Per_correct_mean_D_excl{1,D}(C,:)','descend');
        Top_Feats{C} = I(1:TOP);

        % Get the name and keyword associated with these features
        Top_name(1:TOP) = CodeString(Top_Feats{C},1);
        Top_key(1:TOP) = Keywords(Top_Feats{C},1);
        Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
        Top_mean(1:TOP) = Per_correct_mean_D_excl{1,D}(C,Top_Feats{C});

        Top_40{D}{C} = [Top_name' Top_key' num2cell(Top_ID') num2cell(Top_mean')];

    end   
   
end



%% How many top features of a spe classifier are present across all datasets?

% FIRST, if want to plot more that 40 features, need to run section above
% and change variable TOP

% load Matrix with TopFeat for each classifier and dataset
% load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Top_Feat_Data_Class')  % all specifically removed featured across datasets (combined ('unique'))
whichClassif = 5;
TOP = 500;  % If other than 40, you need to run section above and comment the load (2 lines above)


% Get the top Feat ID of a classifier across datasets
for D = 1:12
    Top_Feat_Cl(:,D) = Top_40{1,D}{1,whichClassif}(:,3);   % 3 refers to 3rd col (Feat ID)
end

Top_Feat_Cl = cell2mat(Top_Feat_Cl);

% Get all elements and assign the number of datasets that share a feature
UniqueElem = unique(Top_Feat_Cl);       
Ncount = histc(Top_Feat_Cl, UniqueElem);   

% Get how much each element is repeated across datasets
for i = 1:size(Ncount,1)
    SumElem(i,1) = sum(Ncount(i,:));
end

% 1st col = ID of top feature, 2nd col = how many datasets have it
Top_Elem = [UniqueElem SumElem];    

% Sort features from most repeated across datasets
[~,I] = sort(SumElem,'descend');
Top_Elem = Top_Elem(I,:);
UniqueElem = UniqueElem(I,:);

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

load('HCTSA.mat', 'Operations')   % Simply to get the list of 7749 features 

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat); % get the equivalent ranking of top feat after special-value features removed
Operations = Operations(equi_Top_Feat,:);  % get Operations of well-behaved features only

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';


% Convert 7749-based index (UniqueElem) to 5603-based index (Idx_UniqueElem)
for i = 1:length(UniqueElem)
    Idx_UniqueElem(i,:) = find(cell2mat(YLabel) == UniqueElem(i));
end

Top_Keywords = Keywords(Idx_UniqueElem(:,1));
Top_name = CodeString(Idx_UniqueElem(:,1));
Top_Keywords = [Top_name Top_Keywords num2cell(Top_Elem)];   % This top features may include features that are top in some datasets and special-value features in other datasets (that's why mean accuracy can be super low (around 7% if 11 datasets produced special values))


% Get % accuracy for each top feat (average of accuracy from ALL datasets)
Accuracy = [];
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')

for ID = 1:length(UniqueElem)
    
    for D = 1:12
        Accuracy(D) = Per_correct_mean_D_excl{1,D}(whichClassif,Idx_UniqueElem(ID));
    end
    
    Top_Keywords(ID,5) = {mean(Accuracy)}; 
    Accuracy = [];

end
    

% Sort from feature generating highest to lowest accuracy across datasets
%%% This is good to know whether some top features, which were considered
%%% 'top' because they yielf the highest accuracy, are actually shared by
%%% most datasets or not. 
[~,Y] = sort(cell2mat(Top_Keywords(:,5)),'descend');
Top_Keywords = Top_Keywords(Y,:);   %  % Features sorted based on their mean accuracy across datasets
Top_ID = Top_Keywords(:,3); 

%%% Most top 3 features are actually present in the top-40-ranking of
%%% approximately 2 to 4 datasets. 


%%% Highlight which top feature was selected from AveragedMatrix

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Get 'Operations of WB features only
load('HCTSA.mat', 'Operations')  

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat);     % remove SV features
Operations = Operations(equi_Top_Feat,:);                 % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(AveragedMatrix_excl(C,:)','descend');
    Top_Feats{C} = I(1:TOP);   
                     
    % Get the name and keyword associated with these features
    Top_name_aver(1:TOP) = CodeString(Top_Feats{C},1);
    Top_key_aver(1:TOP) = Keywords(Top_Feats{C},1);
    Top_ID_aver(1:TOP) = YLabel(Top_Feats{C},1);
    Top_mean_aver(1:TOP) = AveragedMatrix_excl(C,Top_Feats{C});

    Top_40_aver{C} = [Top_name_aver' Top_key_aver' Top_ID_aver' num2cell(Top_mean_aver')];

end


% Top_features of corresponding classifier; used in AverageMatrix
Top_aver = Top_40_aver{1,whichClassif};

% Find how many top features used in AverageMatrix are 
selected_feat_aver = find(ismember(cell2mat(Top_ID(1:TOP)),cell2mat(Top_aver(:,3))) == 1);
% ID_selected = Top_aver(selected_feat_aver,3);

% Plot distribution of top features across datasets
Featmode = cell2mat(Top_Keywords(:,4));   % y axis bar(Featmode)
Features = Top_Keywords(:,1);
Top_ID = Top_Keywords(:,3); 

figure;
h = bar(Featmode);

ax = gca;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:10:length(Features),'uni',0);
xticks(1:10:length(Features))
ylim([0 max(Featmode)+1])
xlabel('Top 500 features across all datasets (ranked from feature yielding highest to lowest accuracy) ')
ylabel('Number of datasets sharing the top feature')
title('Distribution of top features across datasets')

hold on
ii = bar(Featmode(selected_feat_aver));  % Features selected by AveragedMatrix are in red
set(ii,'FaceColor','r');
hold off

% To check why features from AveragedMatrix_excl are not in any top 40
% feature of any dataset
% for D = 1:12
% X(D) = Per_correct_mean_D_excl{1,D}(whichClassif,find(cell2mat(YLabel) == 6660));  % F being a feature in top_40 aver that is not found in UniqueElem
% end
% mean(X)


%%%%% Try something similar, but instead of plotting ind features; let's plot the top families

% Run section above to get Top_Keywords

% Get all Keywords....
Keywords = Top_Keywords(:,2);
UniqueKey = unique(Keywords);
Key_acc = []; 

% And find the index of features belonging to each keyword
for U = 1:length(UniqueKey)
    Key_idx{1,U} = find(strcmp(Keywords,UniqueKey(U)));
    Key_acc = [Key_acc Top_Keywords(Key_idx{1,U},5)];    % The accuracies given by each feature of a family
    Key_idx{2,U} = mean(cell2mat(Key_acc));              % The mean accuracy for each family
    Key_idx{3,U} = length(Key_idx{1,U});
    Key_acc = [];
end

UniqueKey_rep = cell2mat(Key_idx(3,:))';
Mean_acc = cell2mat(Key_idx(2,:))';
UniqueKey = [UniqueKey num2cell(Mean_acc) num2cell(UniqueKey_rep)];   % 3rd col = number of all different features belonging to that family and present in the top 40 features across  datasets

% sort from family having highest to lower mean accuracy (across their
% corresponding features)
[~,I] = sort((Mean_acc)','descend');
UniqueKey = UniqueKey(I,:);


% Plot distribution of top 20 families across datasets
Fammode = cell2mat(UniqueKey(1:20,3));   % y axis bar(Fammode)
Families = UniqueKey(1:20,1)';

figure;
g = bar(Fammode);

ax = gca;
ax.XTickLabels = arrayfun(@(a)string(a),Families(:,1:20),'uni',0);
xticks(1:length(Families))
xtickangle(35)
ylim([0 max(Fammode)+1])
xlabel('Top 20 families across all datasets (ranked from family yielding highest to lowest accuracy) ')
ylabel('Number of datasets sharing one top feature of the family')
title('Distribution of top families across datasets')

yyaxis right
plot(cell2mat(UniqueKey(1:20,2)),'LineWidth',2.5);
ylabel('Accuracy')
%%%% Take example of top feature for given dataset/class and compare it to
%%%% accuracies across other datasets

whichClassif = 1;
% Dataset 1 -> Top feat = [390] IN_AutoMutualInfoStats_diff_20_gaussian.ami14', 90.91 %
% Does feat 390 generate 90% accuracy in all other datasets for Classif 1 ?

for D = 1:12
X(D) = Per_correct_mean_D_excl{1,D}(whichClassif,find(cell2mat(YLabel) == 390));  % F being a feature in top_40 aver that is not found in UniqueElem
end

figure; bar(X)

hold on
ii = bar(X(1));  % Features selected by AveragedMatrix are in red
set(ii,'FaceColor','r');
hold off

xlabel('Datasets'); ylabel('Accuracy'); 
title('Accuracy across datasets for feature 390 (top feature of Dataset 1)')



%% Plot correlation coeff for given classifier across datasets

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')

Combination = {[1 2],[1 3],[1 4],[1 5],[1 6],[1 7],[1 8],[1 9],[1 10],[1 11],[1 12],...
    [2 3],[2 4],[2 5],[2 6],[2 7],[2 8],[2 9],[2 10],[2 11],[2 12],...
    [3 4],[3 5],[3 6],[3 7],[3 8],[3 9],[3 10],[3 11],[3 12],...
    [4 5],[4 6],[4 7],[4 8],[4 9],[4 10],[4 11],[4 12],...
    [5 6],[5 7],[5 8],[5 9],[5 10],[5 11],[5 12],...
    [6 7],[6 8],[6 9],[6 10],[6 11],[6 12],...
    [7 8],[7 9],[7 10],[7 11],[7 12],...
    [8 9],[8 10],[8 11],[8 12],...
    [9 10],[9 11],[9 12],...
    [10 11],[10 12],...
    [11 12]};

for Classif = 1:10   % Choose which classifier to plot
    
    for i = 1:66

        Combi = Combination{i};

        % Get the 5603 accuracies for a given classifier, for both datasets
        Acc_D1_C = Per_correct_mean_D_excl{1,Combi(1)}(Classif,:);  % Accuracy dataset 1
        Acc_D2_C = Per_correct_mean_D_excl{1,Combi(2)}(Classif,:);  % Accuracy dataset 2

        Coefficient = corrcoef(Acc_D1_C,Acc_D2_C);
        Coeff(Classif,i) = Coefficient(2,1);    % Get the correlation coefficient 

    end
end

% If you want to plot corr coeff for all classifiers
Coeff = [Coeff(1,:) Coeff(2,:) Coeff(3,:) Coeff(4,:) Coeff(5,:) Coeff(6,:) Coeff(7,:) Coeff(8,:) Coeff(9,:) Coeff(10,:)];


% Plot correlation coefficients 
figure; histfit(Coeff,10,'kernel')
MEAN = mean(Coeff);
d = xline(MEAN,'k--','mean','LineWidth',2.5);
xlabel('correlation coefficient')
ylabel('distribution of corr coeff between all 12 datasets (66 combinations)')
title('distribution of correlation coefficients')

% Plot scatter plot to visualize correlation for a given classifier
% figure; scatter(Acc_D1_C,Acc_D2_C)


%% Plot correlation coefficients of accuracy between classifiers

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')

Acc_N1N2 = AveragedMatrix_excl(5,:);
Acc_N1REM = AveragedMatrix_excl(7,:);
Acc_N2REM = AveragedMatrix_excl(9,:);

%%% Accuracy N1/N2 vs. N1/REM
Coeff = corrcoef(Acc_N1N2,Acc_N1REM); 
X = Coeff(2,1)   % 0.56

figure; subplot(1,3,1) 
scatter(Acc_N1N2,Acc_N1REM)
xlabel({'N1 vs N2';'61%'}); ylabel({'N1 vs REM                    ';'54%                    '},'Rotation',0)

%%% Accuracy N1/N2 vs. N2/REM
Coeff = corrcoef(Acc_N1N2,Acc_N2REM); 
Y = Coeff(2,1)   % 0.70

subplot(1,3,2) 
scatter(Acc_N1N2,Acc_N2REM)
xlabel({'N1 vs N2';'61%'}); ylabel({'N2 vs REM                    ';'58%                    '},'Rotation',0)

%%% Accuracy N1/REM vs. N2/REM
Coeff = corrcoef(Acc_N1REM,Acc_N2REM); 
Z = Coeff(2,1)  % 0.45

subplot(1,3,3) 
scatter(Acc_N1REM,Acc_N2REM)
xlabel({'N1 vs REM';'54%'}); ylabel({'N2 vs REM                    ';'58%                    '},'Rotation',0)

figure; plot(Acc_N1N2)
hold on; plot (Acc_N1REM)
hold on; plot(Acc_N2REM)
xlim([0 5603])
 
CL = [{Acc_N1N2} {Acc_N1REM} {Acc_N2REM}];
for i = 1:3
    x = CL{i}';
    S = numel(x);
    xx = reshape(x(1:S - mod(S, 20)), 20, []);
    y(:,i)  = sum(xx, 1).' / 20;
end
 
Acc_N1N2 = y(:,1)'; Acc_N1REM = y(:,2)'; Acc_N2REM = y(:,3)'; 
figure; plot(Acc_N1N2,'LineWidth',2)
hold on; plot (Acc_N1REM,'LineWidth',2)
hold on; plot(Acc_N2REM,'LineWidth',2)
xlim([0 265])


%%  Significantly different? 

entropy = [68.8 70.3 68.6 68.4 69.2 66.6 67.9 69.6 66 65.9 67.2 65.8 66];
spectral = [67.2 67.8 67.3 68.7 66.9 66 66.2 66.8];

distr = [90.3 90 89.3 87.3 82.1 81.9 83.3 83.9 84 82.5 84 83.8 82.3];
spectral = [85.2 85.2 86.6 88.6 83.4];

[h,p,ci,stats] = ttest2(entropy,spectral)


for D = 1:12
    M(:,D) = Per_correct_mean_D_excl{1,D}(:,3563);
end
mean(M)

%% Compute pairwise correlation

Top40 = cell2mat(Top_40{1,1}(:,4))';
rho = corr(Top40);

%% Top Features for OVA

% 1) Reconstruct 7749 matrix

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset

InsertCol = zeros(5,1);

Subs = {'001'};

for D = 1:length(Subs)  
    
    sub = Subs{D};
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_OVA(Dataset 001)')

    spec_common = sort([val specifically_removed{1,D}]);
    
    for x = 1:length(spec_common)

        Part1 = Per_correct_mean(:,1:spec_common(x)-1);
        Part2 = [InsertCol Per_correct_mean(:,spec_common(x):end)];
        Per_correct_mean = ([Part1 Part2]);

    end
    
    % Store
    Per_correct_mean_D{D} = Per_correct_mean;
    
end

% 2) Remove all SV features

Per_correct_mean(:,spec_and_common_feat) = [];

% 3) Get top features

% Equivalent indices for WB-features-only data
load('HCTSA.mat', 'Operations')  

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat);   % remove SV features
Operations = Operations(equi_Top_Feat,:);               % get Operations with WB features only

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

% Features reordering: from best to worst feature

Classif = [1:5];

for C = 1:5
    
    Classif = C;   
    
    means = AveragedMatrix(Classif,:);  
    [~,I] = sort((means)','descend');
    AveragedMatrix = AveragedMatrix(:,I);
    
    Top_Feat = I(1:40);   

    Top_name = CodeString(Top_Feat,1);
    Top_key = Keywords(Top_Feat,1);
    Top_ID = YLabel(Top_Feat,1);
    Top_mean = AveragedMatrix(Classif,1:40)';

    Top_40{C} = [Top_name Top_key Top_ID num2cell(Top_mean)];

end




%% Plot accuracy + distribution for AveragedMatrix


% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Get 'Operations of WB features only, from HCTSA_N
load HCTSA.mat   
Operations_ID = [Operations.ID].';

equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); % remove SV features
Operations = Operations(equi_Top_Feat,:);                      % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

TOP = 40;

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(AveragedMatrix_excl(C,:)','descend');
    Top_Feats{C} = I(1:TOP);   
                    
    % Get the name and keyword associated with these features
    Top_name(1:TOP) = CodeString(Top_Feats{C},1);
    Top_key(1:TOP) = Keywords(Top_Feats{C},1);
    Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
    Top_mean(1:TOP) = AveragedMatrix_excl(C,Top_Feats{C});

    Top_40{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];
end


% Indices of classifiers
wake = 0; N1 = 1; N2 = 3; N3 = 3; rem = 5;
CLASSIFIER = {[wake,N1] [wake,N2] [wake,N3] [wake,rem] [N1,N2] [N1,N3] [N1,rem] [N2,N3] [N2,rem] [N3,rem]};
classifstr = {'Wake vs N1','Wake vs N2','Wake vs N3','Wake vs REM','N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};  

%%%
figure;    
[ha, pos] = tight_subplot(1,1,[.065 .05],[.1 .05],[.05 .05]);

% diagonal position og subplots
diag_pos = [1 2 3 4 5 6 7 9 10 13];

% Select a classifier 
for C = 1:1
        
    ClNum = C;  
    Classif = Top_40{1,ClNum};  % Get the top features of the chosen classifier

    
    % Get keywords
    Keyword = char(Classif(:,2));
    
    for k = 1:size(Keyword,1)
        
        F = find(ismember(Keyword(k,:),',') == 1);
        
          if ~isempty(F)
              KEY(k,:) = convertCharsToStrings(Keyword(k,1:F-1));
          else
              KEY(k,:) = convertCharsToStrings(Keyword(k,:));
          end
          
    end
    
        
    % Remove blank spaces
    KEY = strtrim(KEY);

    ALLKEYS{C} = KEY;  % Used for later

    % Get repetitions
    [UniqFam,~,Num] = unique(KEY);
    
    % Compute accuracy for each family
    for i = 1:length(UniqFam)
        Ncount(i) = length(find(UniqFam(i) == KEY));
        idx = find(KEY == UniqFam(i));
        all_accu{i} = Classif(idx,4);
        Accu(i) = mean(cell2mat(all_accu{1,i}));
    end

    % Reorder from highest to lowest accuracy
    [~,I] = sort((Accu)','descend');
    Accu = Accu(:,I);
    Ncount = Ncount(:,I);
    UniqFam = UniqFam(I,:);

    % Later, to know bin edges
    Sum = cumsum(Ncount);
    edges = [0 Sum];

    %%% Plot (rectangle)
    
    [cb] = cbrewer('qual', 'Set3', 12, 'pchip');
    [cb2] = cbrewer('qual', 'Set2', 12, 'pchip');
    [cb] = [cb;cb2];
    
   
    axes(ha(diag_pos(C))); 
    
    for j = 1:length(Accu)
        
        % Color
        load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Top_Features/all_unique_families.mat')
        Color = cb(find(UniqFam(j) == all_unique_families),:);
        
        edges = edges+1;
        pos = [edges(j) 0 Ncount(j) Accu(j)];
        rectangle('Position',pos,'FaceColor',Color)
        text(pos(1)+pos(3)/2,pos(4)+0.8,UniqFam(j),'HorizontalAlignment','left','Rotation',45,'FontSize',12)
        axis([0 60 0 95]) 
    end

    ylim([min(Accu)-4 max(Accu)+4])
    xlim([-1.5 max(edges)+4])
    
    ax=gca;
    yticklabels('auto') 
    
    xlabel('Families of Top 40 features','FontSize',13)
    ylabel('Accuracy (%)','FontSize',13)
    title(classifstr(C),'FontSize',15)
    
    ax.YGrid = 'on';

    clear Accu Ncount KEY
end

delete(ha(8,1))
delete(ha(11,1))
delete(ha(12,1))
delete(ha(14,1))
delete(ha(15,1))
delete(ha(16,1))


%%% Get same colors across datasets
% ALL = [ALLKEYS{1,1};ALLKEYS{1,2};ALLKEYS{1,3};ALLKEYS{1,4};ALLKEYS{1,5};ALLKEYS{1,6};ALLKEYS{1,7};ALLKEYS{1,8};ALLKEYS{1,9};ALLKEYS{1,10}];
% all_unique_families = unique(ALL);

%%% Percentage spectral features among top 40 feat across all classifiers
for x = 1:10
    idx(x) = numel(find(ALLKEYS{1,x} == 'spectral'));
end
PercentSpectral = sum(idx)/(40*10); % 0.20 %



%% Pie chart

%%%% Get top 40 features 

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')


% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Get 'Operations of WB features only, from HCTSA_N
load HCTSA.mat   
Operations_ID = [Operations.ID].';

equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); % remove SV features
Operations = Operations(equi_Top_Feat,:);                      % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

TOP = 40;

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(AveragedMatrix_excl(C,:)','descend');
    Top_Feats{C} = I(1:TOP);   
                    
    % Get the name and keyword associated with these features
    Top_name(1:TOP) = CodeString(Top_Feats{C},1);
    Top_key(1:TOP) = Keywords(Top_Feats{C},1);
    Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
    Top_mean(1:TOP) = AveragedMatrix_excl(C,Top_Feats{C});

    Top_40{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];
end


% Indices of classifiers
wake = 0; N1 = 1; N2 = 3; N3 = 3; rem = 5;
CLASSIFIER = {[wake,N1] [wake,N2] [wake,N3] [wake,rem] [N1,N2] [N1,N3] [N1,rem] [N2,N3] [N2,rem] [N3,rem]};
classifstr = {['Wake vs N1';'          '],['Wake vs N2';'          '],['Wake vs N3';'          '],['Wake vs REM';'           '],['N1 vs N2';'        '],['N1 vs N3';'        '],['N1 vs REM';'         '],['N2 vs N3';'        '],['N2 vs REM';'         '],['N3 vs REM';'         ']};  


figure;    
[ha, pos] = tight_subplot(4,4,[.075 .08],[.1 .05],[.05 .05]);

% diagonal position og subplots
diag_pos = [1 2 3 4 5 6 7 9 10 13];

%%% Plot pie chart: % representation of each family

for C = 1:10

    ClNum = C;  
    Classif = Top_40{1,ClNum};  % Get the top features of the chosen classifier

    % Get keywords
    Keyword = char(Classif(:,2));

    for k = 1:size(Keyword,1)

        F = find(ismember(Keyword(k,:),',') == 1);

        if ~isempty(F)
            KEY(k,:) = convertCharsToStrings(Keyword(k,1:F-1));
        else
            KEY(k,:) = convertCharsToStrings(Keyword(k,:));
        end

    end

    % Remove blank spaces
    KEY = strtrim(KEY);

    % Get repetitions
    [UniqFam,~,Num] = unique(KEY);

    for i = 1:length(UniqFam)
        Ncount(i) = length(find(UniqFam(i) == KEY));
    end

    % Pie chart
    axes(ha(diag_pos(C))); 
    
    labels = UniqFam;
    p = pie(Ncount, labels);
    set(p(2:2:end),'FontSize',12);

    
    % Color
    idx = find(UniqFam == 'spectral');
    
    if idx == 1
        actualidx = 1;
    else
        actualidx = (idx*2)-1;
    end
    
    t = p(actualidx);
    t.FaceColor = 'red';
    title(classifstr(C),'FontSize',13)
        
    clear UniqFam Num KEY Keyword Ncount

end

% Delete unused subplots
delete(ha(8,1))
delete(ha(11,1))
delete(ha(12,1))
delete(ha(14,1))
delete(ha(15,1))
delete(ha(16,1))




%% sup vs unsup, all datasets

%%%%%%%%%%%%%%
%%% Unsup each
%%%%%%%%%%%%%%

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};

    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_each/Per_correct_mean(Dataset %s)',sub))

    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/allfeat_removed')
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

    InsertCol = zeros(size(Per_correct_mean,1),1);

    for x = 1:length(removed_feat_idx{1,D})

        Part1 = Per_correct_mean(:,1:removed_feat_idx{1,D}(x)-1);
        Part2 = [InsertCol Per_correct_mean(:,removed_feat_idx{1,D}(x):end)];
        Per_correct_mean = ([Part1 Part2]);
        removed_feat_idx{1,D} = removed_feat_idx{1,D};

    end

    % Remove both the specific and commonly removed features
    Per_correct_mean(:,spec_and_common_feat) = []; 
    
    % Mean for each classifier
    unsup_each(D,:) = mean(Per_correct_mean');
    
end

%%% Average across datasets
Mean_unsup_each = mean(unsup_each);
std_unsup_each = std(unsup_each);

%%%%%%%%%%%%%%
%%% Unsup all
%%%%%%%%%%%%%%

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};

    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_all/PERC_PER_CLASSIF_10iter_Dataset%s',sub))

    Per_correct_means(D,:) = PERC_PER_CLASSIFIER_10iter;
end

%%% Average across datasets
Mean_unsup_all = mean(Per_correct_means);
std_unsup_all = std(Per_correct_means);


%%%%%%%%%%%%%%
%%% Sup all
%%%%%%%%%%%%%%

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};

    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/sup_all/PERC_PER_CLASSIF_10iter_Dataset%s',sub))

    Per_correct_meanss(D,:) = iteration_svm_testing_accuracy_MEAN;
end

%%% Average across datasets
Mean_sup_all = mean(Per_correct_meanss);
std_sup_all = std(Per_correct_meanss);

% Accuu = Accuu(:,I);
% Per_correct_meanss = Per_correct_meanss(:,I);
% Per_correct_means = Per_correct_means(:,I);
% unsup_each = unsup_each(I,:);

% for i = 1:10
%     [h,p] = ttest2(unsup_each(:,i),Per_correct_means(:,i))
% end

%%%%%%%%%%%%%%
%%% Unsup catch22
%%%%%%%%%%%%%%

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};

    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_catch22/PERC_PER_CLASSIF_10iter_Dataset%s',sub))

    Per_correct_means_c(D,:) = PERC_PER_CLASSIFIER_10iter;
end

%%% Average across datasets
Mean_unsup_catch22 = mean(Per_correct_means_c);


%%%%%%%%%%%%%%
%%% Sup catch22
%%%%%%%%%%%%%%

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};

    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/sup_catch22/PERC_PER_CLASSIF_10iter_Dataset%s',sub))

    Per_correct_means_cc(D,:) = iteration_svm_testing_accuracy_MEAN;
end

%%% Average across datasets
Mean_sup_catch22 = mean(Per_correct_means_cc);


%%%%%%%%%%%%%%
%%% Unsup - Top 5%
%%%%%%%%%%%%%%

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for C = 1:10
    
    for D = 1:length(Subs)  

        sub = Subs{D};

        load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_Top_5perc/average_univariate(280)/Per_correct_mean(Dataset %s)',sub))

        Per_correct_means_five(D,:) = mean(Per_correct_mean(C,:));
    end
    
    %%% Average across datasets
    Mean_uns_Top5perc(C) = mean(Per_correct_means_five);
    

end

%%%%%%%%%%%%%%
%%% Unsup - Top Feat per classifier
%%%%%%%%%%%%%%

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/TopFeatperClass')
Mean_uns_TopFeat = TopFeatperClass(:,2)';

% For std
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ID5603')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')

for i = 1:10
equi_feat(i) = find(TopFeatperClass(i,1)' == ID5603);
end

for C = 1:10   
    for D = 1:12      
        Accuu(D,C) = Per_correct_mean_D_excl{1,D}(C,equi_feat(C));   
    end      
end
    
std_uns_TopFeat = std(Accuu);



% Sort from best to worst classifier (based on unsup all)
[~,I] = sort(Mean_unsup_all,'ascend');
uns_all = Mean_unsup_all(I);
uns_one = Mean_unsup_each(I);
sup_all = Mean_sup_all(I);
uns_catch22 = Mean_unsup_catch22(I);
uns_topFeat = Mean_uns_TopFeat(I);

std_sup_all = std_sup_all(I);
std_uns_TopFeat = std_uns_TopFeat(I);
std_unsup_all = std_unsup_all(I);
std_unsup_each = std_unsup_each(I);

%%%%%%%%%%%%%%%%%%%%%%%
%%% Bar chart, average across pairs
%%%%%%%%%%%%%%%%%%%%%%%

addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';
[GREEN] = cbrewer('div', 'BrBG', 12, 'pchip');
GR_dark = GREEN(12,:);
GR_light = GREEN(10,:);
YELL = GREEN(3,:);
[PURPLE] = cbrewer('div', 'PRGn', 12, 'pchip');
PUR = PURPLE(2,:);
[REDI] =  cbrewer('div', 'RdGy', 12, 'pchip');
RED = REDI(2,:);


Mean_uns_all = mean(uns_all);
Mean_sup_all = mean(sup_all);
Mean_uns_catch22 = mean(uns_catch22);
Mean_uns_one = mean(uns_one);
Mean_uns_TopFeat = mean(uns_topFeat);

% figure; ax=gca;
% 
% BARS = [Mean_sup_all Mean_uns_TopFeat Mean_uns_all Mean_uns_one];
% 
% B = bar(BARS);
% 
% xticklabels({'','','',''})
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% yticklabels(string(0:10:100))
% ylabel('Classification accuracy (%)')
% 
% %%% Just for legend
% hold on;
% a = bar([NaN],'FaceColor',RED); 
% hold on
% b = bar([NaN],'FaceColor',GR_light); 
% hold on;
% c = bar(NaN,'FaceColor',YELL);
% 
% %%% Color bars
% B.FaceColor = 'flat';
% B.CData(1,:) = GR_dark;
% B.CData(2,:) = RED;
% B.CData(3,:) = GR_light;
% B.CData(4,:) = YELL;
% 
% legend([B a b c],'Supervised multivariate classification','Unsupervised clustering using top feature','Unsupervised multivariate clustering','Unsupervised univariate clustering', 'Location','eastoutside','FontSize',24)
% legend boxoff        
% 
% ax.FontSize = 24;
% set(gca,'box','off')
% ax.YGrid = 'on';
% 
% % Add error bars
% STD = [std_sup_all;std_uns_TopFeat;std_unsup_all;std_unsup_each];
% ALL = [sup_all; uns_topFeat; uns_all; uns_one];
% 
% meanSTD = mean(STD');
% meanALL = mean(ALL');
% hold on
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(meanALL);
% for i = 1:1
%     q(i,:) = B(i).XEndPoints;
% end
% er = errorbar(q,meanALL,zeros(size(meanSTD)),meanSTD,'k','linestyle','none');
% 
% set(gcf,'color','white')
% addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
% fpath = '/Users/nico/Documents/HCTSA/Analysis/Clusters';
% %export_fig([fpath filesep '4approaches_3'],'-r 300')

%%%%%%%%%%%%%%%%%%%%%%%
%%% Bar chart,  across pairs
%%%%%%%%%%%%%%%%%%%%%%%
%

ALL = [sup_all; uns_topFeat; uns_all; uns_one];
V = bar(ALL');

figure; ax = gca;

% Add error bar std (plot it first to hide lower bar and lower end
STD = [std_sup_all;std_uns_TopFeat;std_unsup_all;std_unsup_each];

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ALL);
for i = 1:4
    z(i,:) = V(i).XEndPoints;
end
er = errorbar(z,ALL,STD,'k','linestyle','none','LineWidth',1);

hold on
V = bar(ALL');

xlabel('Classifiers')
ylabel('Classification accuracy (%)')

ax = gca;
ax.XTickLabels = {'W/N1','W/N2','N3/W','W/REM','N1/N2','N3/N1','N1/REM','N3/N2','N2/REM','N3/REM'};
ax.XTickLabels = ax.XTickLabels(I);
ax.YTick = 50:5:100;
ax.YTickLabels = string(50:5:100);

xtickangle(30)
ylim([45 101])
ax.FontSize = 12;
set(gca, 'TickLength',[0 0])

%%% Color bars
COLOR = [{GR_dark} {RED} {GR_light} {YELL}];
for B = 1:4
    V(B).FaceColor = 'flat';
    V(B).CData = COLOR{B};
end

%%% Just for legend
hold on
Z = bar([NaN],'FaceColor',GR_dark); 
hold on
a = bar([NaN],'FaceColor',RED); 
hold on
b = bar([NaN],'FaceColor',GR_light); 
hold on
c = bar(NaN,'FaceColor',YELL);

legend([Z a b c],'Supervised feature-based classification','Top-feature univariate clustering','Feature-based clustering','Mean univariate clustering', 'Location','southoutside')
legend boxoff
set(gca,'box','off')
ax.Position = [0.06,0.2,0.9,0.78];
ax.XGrid = 'off';
ax.YGrid = 'on';

set(gca,'fontsize',25)
set(gcf,'color','white')
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'

    fpath = '/Users/nico/Documents/HCTSA/Analysis/Clusters';
   %export_fig([fpath filesep '4approaches_Across_pairs_3'],'-r 300')


%% Mean rank of ind features

% Load data
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% All datasets
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
    
for D = 1:12

    sub = Subs{D};

    % Rank features and get their ID, names, accuracy
    Per_correct_mean = Per_correct_mean_D_excl{1,D};

    means = mean(Per_correct_mean);
    [~,I] = sort((means)','descend');
    Per_correct_mean = Per_correct_mean(:,I);
    Top_Feat = I(1:5603); 

     % Go to corresponding folder and load data
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))      
    load('HCTSA_N.mat', 'Operations')

    % From Operations, take only well behaved features
    equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
    Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
    Operations = Operations(Idx_WB_Feat,:);                      % 'Operations' with WB features only

    % Get the name and keyword associated with these features
    CodeString = Operations.CodeString;  
    Keywords = Operations.Keywords;
    YLabel = Operations.ID;

    Top_name(1:5603) = CodeString(Top_Feat,1);
    Top_key(1:5603) = Keywords(Top_Feat,1);
    Top_ID(1:5603) = YLabel(Top_Feat,1);

    Top_5603 = [Top_name' Top_key' num2cell(Top_ID')];

    % Top_mean (mean over classifiers)
    for F = 1:length(Top_5603)
        for C = 1:10
            Accu(F,C) = mean(Per_correct_mean(C,F));
        end
    end

    Top_mean = mean(Accu');
    Top_5603 = [Top_5603 num2cell(Top_mean')];
    
    Rank(:,D) = Top_5603(:,3);   % Rank of all features for all datasets
    
    clear top_5603
    
end


Rank = cell2mat(Rank);

% Now, get the mean rank for each feature
Feat = sort(Top_ID, 'ascend');

for F = 1:5603
    
    for D = 1:12
    
        FeatID(D) = find(Feat(F) == Rank(:,D));
    
    end
    
    % Mean rank across datasets
    ID_allD(F) = mean(FeatID);
    
    FeatID = []; % clear var
    
end

% Sort from lowest to highest rank (consistently best to worst)
[~,I] = sort((ID_allD)','ascend');   % I = top rank (/5603)
ID_allD_rank = ID_allD(I);

% Convert rank (/5603  ->  /7749)
Real_ID = Feat(I);

% Now get ID, name and keyword
Top_10_ID = Real_ID(1:10);
Top_10_KEY = Keywords(I(1:10));
Top_10_NAME = CodeString(I(1:10));

Top_10 = [num2cell(Top_10_ID)' Top_10_KEY Top_10_NAME num2cell(ID_allD_rank(1:10))'];

%%% Get the accuracy as well
for F = 1:10
    for D = 1:12
        Accuracy(D,F) = mean(Per_correct_mean_D_excl{1,D}(:,I(F)));
    end
end

Mean_Accu_acr_D = mean(Accuracy);
%STD across datasets for each feature
Accu_STD = std(Accuracy);

%%% Get ranks across datasets for top 10 features
for F = 1:10
    for D = 1:12
        ranks(D,F) = find(Top_10_ID(F) == Rank(:,D));
    end
end

%%%%%% Get z_score 
for D = 1:12     % For each dataset. Take one feature, get the z score (perf relative to all other features)
    
    ACCU = mean(Per_correct_mean_D_excl{1,D});   % Take accuracy (averaged across classifiers)
    Z(D,:)= zscore(ACCU); 
    
end

% Sort from highest to lowest Z
Mean_Z = mean(Z);
 [~,I] = sort(Mean_Z,'descend');
Mean_Z = Mean_Z(I);   

% Get their names and stuff
TOP = 5603;

Best10 = I(1:TOP);
Top_10_ID = Operations.ID(Best10);
Top_10_KEY = Operations.Keywords(Best10);
Top_10_NAME = Operations.CodeString(Best10);

%%% Get the accuracy as well
for F = 1:TOP
    for D = 1:12
        Accuracy(D,F) = mean(Per_correct_mean_D_excl{1,D}(:,I(F)));
    end
end

Mean_Accu_acr_D = mean(Accuracy);
Accu_STD = std(Accuracy);
Top_Mean = Mean_Accu_acr_D(1:TOP);

Top_10 = [num2cell(Top_10_ID) Top_10_KEY Top_10_NAME num2cell(Top_Mean)' num2cell(Mean_Z(1:TOP))' ];

% Get Best 5%
Top_5perc = Top_10;


%%% Plot boxplot
figure; boxplot(ranks)
title('Median rank of top 10 features')
xlabel('Features')
ylabel('Median rank')
    set(gcf,'color','white')
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
%     export_fig([fpath filesep 'BoxplotTop10'],'-r 300')

%%% Std 
for i = 1:10
    STD(i) = std(ranks(:,i));
end

%%% Plot the std

% Plot std (variance across classifiers)
figure; ax = gca;
stdshade(Accuracy,0.5,'r');
title('Mean accuracy across features')
xlabel('features sorted from best to worst')
ylabel('Percentage Accuracy')


%% Mean rank of best feature for each classifier (1800005)

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'005'};  % '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
end

Per_correct_mean = Per_correct_mean_D_excl{y};

% Load Data
load HCTSA_N.mat

% From Operations, take only well behaved features
equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
Operations = Operations(Idx_WB_Feat,:);  

% Get the name and keyword associated with these features
CodeString = Operations.CodeString;  
Keywords = Operations.Keywords;
YLabel = Operations.ID;

TOP = 5603;

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(Per_correct_mean(C,:)','descend');
    Top_Feats{C} = I(1:TOP);   
                     
    % Get the name and keyword associated with these features
    Top_name(1:TOP) = CodeString(Top_Feats{C},1);
    Top_key(1:TOP) = Keywords(Top_Feats{C},1);
    Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
    Top_mean(1:TOP) = Per_correct_mean(C,Top_Feats{C});

    Top_5603{C} = [Top_name' Top_key' num2cell(Top_ID') num2cell(Top_mean')];

end


% Index of classifier involving each stage
Wake = [1 2 3 4];    
N1 = [1 5 6 7];    
N2 = [2 5 8 9];    
N3 = [3 6 8 10];    
REM = [4 7 9 10];    

ST = [{Wake} {N1} {N2} {N3} {REM}];

% Now for each stage, get the mean rank across the4 corresponding
% classifiers

% load index of 5603 features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ID5603')
 
for F = 1:5603   % for each feature
    
    for S = 1:5    % for each stage 
    
        Classif = ST{S};
        
        for C = 1:4     % for each of the 4 classifier
                        
            Rank = Top_5603{1,Classif(C)}(:,3);
            
            FeatID(C) = find(ID5603(F) == cell2mat(Rank));
            
        end

        % Mean across 4 classifiers (1 stage) for 1 Feature
        Mean_stage(F,S) = mean(FeatID);
        
    end
    
end

for S = 1:5
    
    % Sort from lowest to highest rank (consistently best to worst)
    % Say, N1
    [~,Y] = sort(Mean_stage(:,S)','ascend');   % I = top rank (/5603)
    Mean_stage(:,S) = Mean_stage(Y,S);

    % Now get ID, name and keyword
    Top_10_ID = table2cell(Operations(Y(1:10),4));
    Top_10_KEY = table2cell(Operations(Y(1:10),3));
    Top_10_NAME = table2cell(Operations(Y(1:10),1));

    Top_10 = [Top_10_ID Top_10_KEY Top_10_NAME];

    %%% Get the accuracy as well
    for F = 1:10   % for each feature


        Classif = ST{S};

        for C = 1:4     % for each of the 4 classifier

            IDX = find(cell2mat(Top_5603{1,Classif(C)}(:,3)) == cell2mat(Top_10_ID(F))) ;
            MEANS(C) = Top_5603{1,Classif(C)}(IDX,4);

        end

        % Mean across 4 classifiers (1 stage) for 1 Feature
        THEMEAN_stage(F,1) = mean(cell2mat(MEANS));

    end

Top_10 = [Top_10_ID Top_10_KEY Top_10_NAME num2cell(THEMEAN_stage) num2cell(Mean_stage(1:10,S))];

TOPP{S} = Top_10;

clear Top_10 Top_10_ID Top_10_KEY Top_10_NAME THEMEAN_stage

end

%% Catch 22

% Find the right index of 22 features

FEAT = {'DN_HistogramMode_5','DN_HistogramMode_10','first_1e_ac', ...
    'first_min_acf', 'CO_HistogramAMI_even_5_2','CO_trev_1_num', ...
    'MD_hrv_classic_pnn40','SB_BinaryStats_mean_longstretch1', ...
    'SB_TransitionMatrix_3ac_sumdiagcov','PD_PeriodicityWang_th0.01', ...
    'CO_Embed2_Dist_tau_d_expfit_meandiff','IN_AutoMutualInfoStats_40_gaussian_fmmi', ...
    'FC_LocalSimple_mean1_tauresrat','DN_OutlierInclude_p_001_mdrmd', ...
    'DN_OutlierInclude_n_001_mdrmd','SP_Summaries_welch_rect_area_5_1', ...
    'SB_BinaryStats_diff_longstretch0','SB_MotifThree_quantile_hh', ...
    'SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1','SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1', ...
    'SP_Summaries_welch_rect_centroid','FC_LocalSimple_mean3_stderr'};
                
% Load Operation Names, create separate variable
load('HCTSA.mat','Operations')
Name = {Operations.Name}.';

% Get index (/7749) of catch22 
for F = 1:22
    Index(F) = find(contains(Name,FEAT(F)));
end

%%% 18th feature has 3 'duplicate' (because "contains" gets 3) but its 3604

Catch22_idx = [11,12,134,135,243,1121,7634,3477,1406,1585,1965,310,2997,3264,3294,4492,3467,3604,4036,4156,4421,3010];


%% TRYYY

load('/Users/nico/Documents/HCTSA/Analysis/clustering_using_best_feat/BestFeat_unsup')
load('/Users/nico/Documents/HCTSA/Analysis/clustering_using_best_feat/AllofThem_2')

[~,I] = sort(AllofThem_2,'ascend');
AllofThem_2 = AllofThem_2(I);
PERC_PER_CLASSIFIER_10iter = PERC_PER_CLASSIFIER_10iter(I);

% Line plot
figure; 
h = plot(1:10,AllofThem_2,'LineWidth',1.3,'Color',[0.3010 0.7450 0.9330]);  % light blue
hold on
i = plot(1:10,PERC_PER_CLASSIFIER_10iter,'LineWidth',1.3,'Color','g');  % green 
hold off



% legend('Unsup - using all features','SVM - using all features','Unsup - one feature at a time','SVM - one feature at a time','Location','eastoutside')
xlabel('Classifiers')
ylabel('Classification accuracy (%)')

ax = gca;
ax.XTick = 1:11;
ax.XTickLabels = {'W vs N1','W vs N2','N3 vs W','W vs REM','N1 vs N2','N3 vs N1','N1 vs REM','N3 vs N2','N2 REM','N3 vs REM'};
ax.XTickLabels = ax.XTickLabels(I);
xtickangle(30)
ylim([45 100])
yline(50,'--','chance level','FontSize',18);
ax.FontSize = 12;

legend('Full set of features','Best features','Location','eastoutside')
legend boxoff
set(gca,'fontsize',20)
    set(gcf,'color','white')
%     fpath = '/Users/nico/Documents/HCTSA/Analysis';
%     export_fig([fpath filesep 'sup unsup catch22_2'],'-r 300')


%% Corr Mat and Violin Plot from top 10 Feat across datasets, in one dataset


%%%% CorrMat

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 

    % Load Data
    load('HCTSA_N.mat','Operations')

    % Top 10 Feat across datasets
    Top10Feat = [4449 4325 922 923 2117 2133 921 4444 4320 4343]';
    NAME = {'SPE1' 'SPE2','ENT1','ENT2','STA1','STA2','ENT3','SPE3','SPE4','SPE5'};
    
    % Get equivalent index for corresponding dataset
    for i = 1:10    
        equi_feat(i) = find(ismember(Operations.ID,Top10Feat(i)));
    end

    Top_name = Operations.CodeString(equi_feat);
    Top_key = Operations.Keywords(equi_feat);
    Top_ID = Operations.ID(equi_feat);

    Top_10 = [Top_name Top_key num2cell(Top_ID)];

    % Top_mean 
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_each/Per_correct_mean(Dataset %s)',sub))
    for F = 1:10
        Top_10_Acc(F) = mean(Per_correct_mean(:,equi_feat(F)));
    end

    Top_10 = [Top_10 num2cell(Top_10_Acc')];

end

load('HCTSA_N.mat','TS_DataMat')

% Parameters
numTopFeatures = 10;   
op_ind = equi_feat';      
distanceMetric = 'abscorr';
clusterThreshold = 0.2;

% Compute correlations based on hctsa responses
EEGonly = 1:size(TS_DataMat,1)/7;
TS_DataMat = TS_DataMat(EEGonly,op_ind);    
Dij = BF_pdist(TS_DataMat','abscorr');  

% Ylabels
YLabel = [];
for i = 1:length(op_ind)
    YLabel = [YLabel {sprintf('%s (%1.1f%%)',NAME{i},string(Top_10_Acc(:,i)))}];
end

% Plot the correlation matrix
[~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',YLabel);
set(gca,'fontsize',22)

set(gcf,'color','white')
fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
% export_fig([fpath filesep 'CorrMattop10_on_D1_Paper'],'-r 300')
    
%%%%% Violin plot

ifeat = equi_feat';

load HCTSA_N.mat
% Make data structure for TS_SingleFeature
data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
                'Operations',Operations);
            
% Set the colors to be assigned to groups:
numClasses = 5;
colors = GiveMeColors(numClasses);

% Space the figures out properly:
numClasses = 5; 
subPerFig = 10; % subplots per figure
numFeaturesDistr = 10;  % How many top features I want to plot
numFigs = ceil(numFeaturesDistr/subPerFig);

    % rearrange according to corrmat
    % arr = [1 10 2 3 9 6 7 4 5 8];
    arr = [7 3 4 6 5 9 8 1 2 10];
    equi_feat = equi_feat(arr);
    Top_10 = Top_10(arr,:);
    NAME = NAME(arr);

    ifeat = equi_feat;

    for figi = 1:10

        % Get the indices of features to plot
        r = ((figi-1)*subPerFig+1:figi*subPerFig);

        if figi==numFigs % filter down for last one
            r = r(r<=numFeaturesDistr);
        end

        featHere = ifeat(r); % features to plot on this figure

        % Make the figure
        f = figure('color','w');
        f.Position(3:4) = [1353, 857];
        % Loop through features
        for opi = 1:length(featHere)
            subplot(ceil(length(featHere)/5),5,opi);
            TS_SingleFeature_1D_2(data,featHere(opi),true,false,opi);
           if numel(Top_10{opi,1}) > 25
               Top_10{opi,1} = Top_10{opi,1}(1:25);
           end
            title({sprintf('%s (%1.1f%%)',NAME{opi},string(Top_10_Acc(:,arr(opi))))},'interpreter','none');
           set(gca,'fontsize',16)

        end


    end

        set(gcf,'color','white')
        fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
        export_fig([fpath filesep 'Violintop10_on_D1_Paper'],'-r 300')

%%%%% Accuracy across classifiers
load HCTSA_N.mat
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')

equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
Operations = Operations(Idx_WB_Feat,:);  

CodeString = Operations.CodeString;  
Keywords = Operations.Keywords;
YLabel = Operations.ID;

EQUI = find(ismember(YLabel,Top_ID));
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
for F = 1:10
    for D = 1:12
    FEAT(:,D) = Per_correct_mean_D_excl{1,D}(:,EQUI(F));
    end
    
    Mean_FEAT(F,:) = mean(FEAT');
    
end
    
    
    


%% Violin plot, best feature from each classifier across datasets, on Dataset 1

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
end

Per_correct_mean = Per_correct_mean_D_excl{y};

% Load Data
load HCTSA_N.mat

% From Operations, take only well behaved features
equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
Operations = Operations(Idx_WB_Feat,:);  

% Get the name and keyword associated with these features
CodeString = Operations.CodeString;  
Keywords = Operations.Keywords;
YLabel = Operations.ID;
TOP = 5603;

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:12
    
    sub = Subs{D};
   
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))      

    for C = 1:10

        % Features reordering: from best to worst feature
        [~,I] = sort(Per_correct_mean_D_excl{1,D}(C,:)','descend');
        Top_Feats{C} = I(1:TOP);

        % Get the name and keyword associated with these features
        Top_name(1:TOP) = CodeString(Top_Feats{C},1);
        Top_key(1:TOP) = Keywords(Top_Feats{C},1);
        Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
        Top_mean(1:TOP) = Per_correct_mean_D_excl{1,D}(C,Top_Feats{C});

        Top_40{D}{C} = [Top_name' Top_key' num2cell(Top_ID') num2cell(Top_mean')];

    end   
   
end

%%% For each dataset, each classifier, get rank of each feature
Feat = sort(Top_ID, 'ascend');

%
%%% Get z_score 
TOP = 5603;
    
for C = 1:10
    for D = 1:12     % For each dataset. Take one feature, get the z score (perf relative to all other features)

        ACCU = Per_correct_mean_D_excl{1,D}(C,:);   % Take accuracy (averaged across classifiers)
        Z(D,:) = zscore(ACCU); 

    end
    
    Mean_Z(C,:) = mean(Z);
    
    [~,I] = sort(Mean_Z(C,:),'descend');
    Mean_Z(C,:) = Mean_Z(C,I);   
    
    Best1 = I(1:TOP);
    Top_IDS(C,:) = Operations.ID(Best1);
    Top_KEY(C,:) = Operations.Keywords(Best1);
    Top_NAME(C,:) = Operations.CodeString(Best1);
    
    for y = 1:TOP
        for D = 1:12
           ACCUS(D,y) = Per_correct_mean_D_excl{1,D}(C, Best1(y));
        end
    end
    
    Top_ACCU(C,:) = mean(ACCUS);
    
end

for i = 1:10
    X(i) = find(Top_IDS(i,:) == 4449)
    Y(i) = Top_ACCU(i,X(i))
    
end


TopFeatperClass = [Top_ID(:,1) Top_ACCU];

Classif = {'Wake/N1','Wake/N2','Wake/N3','Wake/REM','N1/N2','N1/N3','N1/REM','N2/N3','N2/REM','N3/REM'};
Table2 = array2table([Classif' num2cell(Top_ID(:,1)) Top_KEY Top_NAME num2cell(Top_ACCU)],'VariableNames',{'Classifier','ID','Keyword','Name','Accuracy'});


%%%%% Violin plot
sub = '001';
cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';
FeatName = {'F1','F2','F3','F4','F5','F6','F7','F8','F9','F10'};


load HCTSA_N.mat
% Make data structure for TS_SingleFeature
data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
                'Operations',Operations);

ID = Operations.ID;
for i =1:10
ifeat(i) = find(ID == Top_ID(i,1));
end

% Set the colors to be assigned to groups:
numClasses = 5;
colors = GiveMeColors(numClasses);

% Space the figures out properly:
numClasses = 5; 
subPerFig = 10; % subplots per figure
numFeaturesDistr = 10;  % How many top features I want to plot
numFigs = ceil(numFeaturesDistr/subPerFig);

for figi = 1:10

    % Get the indices of features to plot
    r = ((figi-1)*subPerFig+1:figi*subPerFig);

    if figi==numFigs % filter down for last one
        r = r(r<=numFeaturesDistr);
    end

    featHere = ifeat(r); % features to plot on this figure

    % Make the figure
    f = figure('color','w');
    f.Position(3:4) = [1353, 857];
    % Loop through features
    for opi = 1:length(featHere)
        subplot(ceil(length(featHere)/5),5,opi);
        TS_SingleFeature_1D_2(data,featHere(opi),true,false,opi);
       if numel(Top_NAME{opi,1}) > 25
           Top_NAME{opi,1} = Top_NAME{opi,1}(1:20);
       end
        title({[Classif{opi}];sprintf('%s',FeatName{opi})},'interpreter','none');
       set(gca,'fontsize',16)

    end


end

set(gcf,'color','white')
fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
%export_fig([fpath filesep 'Violintop1_on_D1_EachClass_Paper'],'-r 300')


%% Only 1 violin plot (4449)
sub = '001'
DD = str2num(sub);

cd('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_001') 
load HCTSA.mat

% Make data structure for TS_SingleFeature
data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
                'Operations',Operations);
            
numClasses = 5; 
subPerFig = 9; % subplots per figure
numFeaturesDistr = 9;  % How many top features I want to plot
numFigs = ceil(numFeaturesDistr/subPerFig);

ifeat = [7160 2780 4449 4343 4008 2762 363 6594 7568];

Top_MEAN = [76.6 76.5 83.2 77.6 76.2 77.1 77.5 77.4 76.7];
for i = 1:length(ifeat)
    Top_NAME(i) = cellstr(Operations(ifeat(i)).CodeString);
    Top_KEY(i) = cellstr(Operations(ifeat(i)).Keywords);
    Top_ID(i) = ifeat(i);
end


for figi = 1:9

    % Get the indices of features to plot
    r = ((figi-1)*subPerFig+1:figi*subPerFig);

    if figi==numFigs % filter down for last one
        r = r(r<=numFeaturesDistr);
    end

    featHere = ifeat(r); % features to plot on this figure

    % Make the figure
    f = figure('color','w');
    f.Position(3:4) = [1353, 857];
    % Loop through features
     for opi = 1:length(featHere)
            subplot(ceil(length(featHere)/5),5,opi);
            TS_SingleFeature_1D_2(data,featHere(opi),true,false,opi,DD);
           if numel(Top_NAME{opi}) > 25
               Top_NAME{opi} = Top_NAME{opi}(1:20);
           end
            title({sprintf('%s',Top_NAME{opi})},'interpreter','none');
        
         set(gca,'fontsize',16)
         set(0,'DefaultAxesTitleFontWeight','normal');

    end

end

f.Position = [395,208,430,589];

addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
set(gcf,'color','white')
fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
% export_fig([fpath filesep 'Top40Feat_violin'],'-r 300')


%% Paper figure: Top 40 features CoorMat and violin plots on D001)

% Load index 5603 features (common to all datasets)
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ID5603')

%%%% CorrMat
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))

    % Load Data
    load('HCTSA_N.mat','Operations')
    
    % Convert number feat (eg 6006) to 5603 WB features
    X = find(ismember(Operations.ID,ID5603));
    
    % Load Per_correct_mean and turn to 5603
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_each/Per_correct_mean(Dataset %s)',sub))
    PER_correct_mean_allD(D,:) = mean(Per_correct_mean(:,X));
    
end
    
% Average across datasets (and classifiers already)
Per_correct_mean = [];
Per_correct_mean = mean(PER_correct_mean_allD);

% Sort from best to worst feature
[~,I] = sort(Per_correct_mean','descend');
Per_correct_mean = Per_correct_mean(I);

% Which dataset are we using to use for feature values?

AllD = {'001','005','439','458','596','748','749','752','604','807','821','870'};
cd('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_001')

load('HCTSA.mat','Operations')

% Get equivalent index for corresponding dataset (whatever the dataset)
ID = [Operations.ID].';
X = find(ismember(ID,ID5603));
Operations = Operations(X,:);

% Get Name, Keyword, ID, Accuracy of top features
TOP = 40;
TopFeat = I(1:TOP);

Top_ID = {Operations(TopFeat,:).ID}.';
Top_Name = {Operations(TopFeat,:).CodeString}.';
Top_Key = {Operations(TopFeat,:).Keywords}.';
Top_Mean = Per_correct_mean(1:TOP)';

Top_Feat = [Top_ID Top_Name Top_Key num2cell(Top_Mean)];

for D = 3
    
Data = AllD{D};
cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',Data))

% Parameters for CorrMat
distanceMetric = 'abscorr';
clusterThreshold = 0.2;

% Load feature values (TS_Datamat) for one dataset
% load HCTSA_N for normalised or HCTSA for unnorm) 
% Using Spearman here (in BF_pdist, change line 123 and replace 'corr 'Pearson' by 'spe'. 
% If use Spearman: no diff between unnorm and norm matrix. If Pearson,
% yeilds different unnorm matrix
load('HCTSA_N.mat', 'Operations','TS_DataMat')
ID2 = [Operations.ID].';
X = find(ismember(ID2,ID5603));
EEGonly = 1:size(TS_DataMat,1)/7;

TS_DataMat = TS_DataMat(EEGonly,X);  % EEG only, with the 5603 WB features  
TS_DataMat = TS_DataMat(:,TopFeat);  % Only the top 40 features

% Compute correlations based on hctsa responses  
Dij = BF_pdist(TS_DataMat','abscorr');  

% Ylabels
YLabel = [];
ZLabel = [];

for i = 1:TOP
    YLabel = [YLabel {sprintf('%s (%1.1f%%)',Top_Name{i}(4:end),string(Top_Mean(i,:)))}];
    ZLabel = [ZLabel {sprintf('%s',Top_Name{i}(1:2))}];
end
 
% Plot the correlation matrix
[~,cluster_Groupi,ord,~,WhichREP,objectLabels] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',YLabel);

%%% Add text with keywords
KEYS = ZLabel(ord);

%%% Colors of text
addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';

l_green = cbrewer('div', 'PiYG', 12, 'pchip'); l_green = l_green(9,:);
d_red = cbrewer('div', 'RdBu', 12, 'pchip'); d_red = d_red(2,:);
grey = cbrewer('div', 'RdGy', 12, 'pchip'); grey = grey(10,:);
l_blue = cbrewer('div', 'RdBu', 12, 'pchip'); l_blue = l_blue(9,:);
brown = cbrewer('div', 'BrBG', 12, 'pchip'); brown = brown(2,:);
purple = cbrewer('div', 'PRGn', 12, 'pchip'); purple = purple(2,:);
d_green = cbrewer('div', 'PRGn', 12, 'pchip'); d_green = d_green(11,:);
d_blue = cbrewer('div', 'RdBu', 12, 'pchip'); d_blue = d_blue(11,:);
orange = cbrewer('div', 'PuOr', 12, 'pchip'); orange = orange(4,:);
pink = cbrewer('div', 'PiYG', 12, 'pchip'); pink = pink(4,:);
yellow =  cbrewer('qual', 'Set3', 12, 'pchip'); yellow = yellow(12,:);

Colors = [{l_green} {d_red} {grey} {d_blue} {brown} {purple} {d_green} {l_blue} {orange} {pink} {yellow}];
All_KEYS = {'SC','SP','SY','NW','IN','WL','EN','FC','PH','NL','MF'};

%%% Plot text
gap = 1;
for i = 1:TOP
    
    index = find(string(KEYS(:,i)) == All_KEYS);
    
    text(-1.4,gap,string(KEYS(:,i)),'FontSize',14,'FontWeight','bold','color',Colors{index})
    gap = gap+1;
end


set(gcf,'color','white')
fpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Top_Features/all_classifiers(ind_dataset)';
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
%export_fig([fpath filesep sprintf('CM_top40(%s)',Data)],'-r 300')

end

%%%%% Violin plot

% remove bar
for i = 1:numel(WhichREP)
WhichREP(1,i) = erase(WhichREP(1,i),'\');
end
% Reorder according to corrMat
WhichREP = WhichREP([5 7 4 6 2 1 3 8]);


for i = 1:numel(WhichREP)
    
    A = find(contains(string(YLabel),string(WhichREP(i))));
    ifeat(i) = Top_ID{A};

end

for i = 1:numel(ifeat)
    V = find(ismember(cell2mat(Top_ID),(ifeat(i))));
    Violin_name(i) = Top_Name(V);
    Violin_accu(i) = Top_Mean(V);
    LenghtName(i) = numel(Violin_name{1,i});
end

load HCTSA.mat  % Switch between HCTSA.mat and HCTSA_N.mat to get (un)norm

% Make data structure for TS_SingleFeature
data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
                'Operations',Operations);
            
% Set the colors to be assigned to groups:
numClasses = 5;
colors = GiveMeColors(numClasses);

% Space the figures out properly:
numClasses = 5; 
subPerFig = 8; % subplots per figure
numFeaturesDistr = 8;  % How many top features I want to plot
numFigs = ceil(numFeaturesDistr/subPerFig);
DD = str2num(Data);

%
for figi = 1:numel(ifeat)

    % Get the indices of features to plot
    r = ((figi-1)*subPerFig+1:figi*subPerFig);

    if figi==numFigs % filter down for last one
        r = r(r<=numFeaturesDistr);
    end

    featHere = ifeat(r); % features to plot on this figure
    ID1 = [Operations.ID].';

    for i = 1:numel(ifeat)
        featHere_D(i) = find(ismember(ID1,featHere(i)));
    end

    % Make the figure
    f = figure('color','w');
    f.Position = [102,52,1120,745];
    
    % Parameters to decide gap between subplots
    % Default first subplot :  0.13 0.5838 0.1566 0.3412
    % Next plots of top row: 0.3361 _ _ _ (+0.2061)
    % Going row down: _ 0.5838 _ _ to _ 0.11 _ _ (-0.4738)


    % Loop through features
    for opi = 1:length(featHere)
        P = subplot(ceil(length(featHere)/4),4,opi);
            
        if  opi ==  2 || opi == 3 || opi == 4 
            P.Position(1) = 0.13+(0.2*(opi-1));
        end
        if opi ==  6 || opi == 7 || opi == 8
           P.Position(1) = 0.13+(0.2*(opi-5));
        end
        if  opi == 5 || opi ==  6 || opi == 7 || opi == 8
            P.Position(2) = 0.13;
        end
        
        TS_SingleFeature_1D_2(data,featHere_D(opi),true,false,opi,DD);
        
       if numel(Violin_name{1,opi}) > 21
          Violin_names{1,opi} = Violin_name{1,opi}(1:21);
          Violin_names{2,opi} = Violin_name{1,opi}(22:end);
       else
          Violin_names{1,opi} = Violin_name{1,opi}(1:end);
       end
       
       T = title({[sprintf('%s',Violin_names{1,opi})];[[sprintf('%s',Violin_names{2,opi})],[sprintf(' (%1.1f%%)',string(Violin_accu(:,opi)))]]},'interpreter','none');
       
       set(gca,'fontsize',16)
       T.FontSize = 15;
       T.FontWeight = 'bold';

    end


end

set(gcf,'color','white')
fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
%export_fig([fpath filesep 'Violintop40_unnorm_439'],'-r 300')





