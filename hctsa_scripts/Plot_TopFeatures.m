

%% Identify the top features

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))
    
    %%% Features reordering;
    means = mean(Per_correct_mean);
    [~,I] = sort((means)','descend');
    Per_correct_mean = Per_correct_mean(:,I);

    %%% Best X features 
    Top_Feat(:,D) = I(1:20);
    
end

gpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Top_50_feat';
save(fullfile(gpath,'Top_Feat.mat'))

%% Features reordering based on TS_CLUSTER (op_clust)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D')

TS_Cluster
Per_correct_mean = Per_correct_mean_D{1,1};
Per_correct_mean = Per_correct_mean(:,(op_clust.ord)');


%% Return the top features from averaged matrix (WB feat)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

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

Top_Feat = I(1:10);   

Top_name = CodeString(Top_Feat,1);
Top_key = Keywords(Top_Feat,1);
Top_ID = YLabel(Top_Feat,1);
Top_mean = mean(AveragedMatrix(:,1:10))';

Top_10 = [Top_name Top_key Top_ID num2cell(Top_mean)];


%% Same as above but top features for each classifier

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

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


%% Get top features EXCLUDING the special values features

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

load('HCTSA.mat', 'Operations')   % Simply to get the list of 7749 features 

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat); % get the equivalent ranking of top feat after special-value features removed
Operations = Operations(equi_Top_Feat,:);  % get Operations of well-behaved features only

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

%%% and replace AverageMatrix and Per_correct_mean_D by AverageMatrix_excl and Per_correct_mean_D_excl


%% Line plot of accuracy from all features sorted from best to worst

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_incl_all_feat_removed(7749)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')  % all specifically removed featured across datasets (combined ('unique'))
AveragedMatrix(:,spec_and_common_feat) = [];

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
figure; plot(MeanFeat,'LineWidth',1.3)
title('Mean accuracy across features')
xlabel('features (sorted from "best" to "worst")')
ylabel('Percentage Accuracy')

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

%% on obtaining Per_correct_mean_D  /   or Per_correct_mean_D_SVM

% Copy pasted script to get 7749-matrix
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset

InsertCol = zeros(10,1);

Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once
% (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering

    % I'll insert both the commonly and specifically removed features for
    % this dataset
    spec_common = sort([val specifically_removed{1,y(D)}]);
    
    for x = 1:length(spec_common)

        Part1 = Per_correct_mean(:,1:spec_common(x)-1);
        Part2 = [InsertCol Per_correct_mean(:,spec_common(x):end)];
        Per_correct_mean = ([Part1 Part2]);

    end
    
    % Store
    Per_correct_mean_D_SVM{D} = Per_correct_mean;
    
end

%% SVM 

%%% 1) Section above ("on obtaining Per_correct_mean_D / or
%%% Per_correct_mean_D_SVM") to get Per_correct_mean_D_SVM
%%% 2) Remove SV features:
%%% Per_correct_mean_D_SVM{1,1}(:,spec_and_common_feat) = []; to get
%%% 10x5308 matrix


%% Top 40 features for one Dataset (only WB Features)

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
end

Per_correct_mean = Per_correct_mean_D_SVM{y};

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

%%%%% Get the AveragedMatrix with all 12 datasets and their 5308 WB features

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

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
HCTSA_N = false;   % true if want to plot normalized hctsa values; else plot HCTSA.mat values

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
title(sprintf('Dependencies between %u top features (organized into %u clusters)',...
                        numTopFeatures,length(cluster_Groupi)))
                    

%% Corr matrix on specific classifiers

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
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

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(Per_correct_mean(C,:)','descend');
    Top_Feats{C} = I(1:40);   
                     
    % Get the name and keyword associated with these features
    Top_name(1:40) = CodeString(Top_Feats{C},1);
    Top_key(1:40) = Keywords(Top_Feats{C},1);
    Top_ID(1:40) = YLabel(Top_Feats{C},1);
    Top_mean(1:40) = Per_correct_mean(C,Top_Feats{C});

    Top_40{C} = [Top_name' Top_key' num2cell(Top_ID') num2cell(Top_mean')];

end


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
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

% Get 'Operations of WB features only, from HCTSA_N
load HCTSA.mat   
Operations_ID = [Operations.ID].';

equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); % remove SV features
Operations = Operations(equi_Top_Feat,:);                      % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

for C = 1:10
    
    % Reorder features: from yielding theg highest to lowest accuracy
    [~,I] = sort(AveragedMatrix_excl(C,:)','descend');
    Top_Feats{C} = I(1:40);   
                     
    % Get the name and keyword associated with these features
    Top_name(1:40) = CodeString(Top_Feats{C},1);
    Top_key(1:40) = Keywords(Top_Feats{C},1);
    Top_ID(1:40) = YLabel(Top_Feats{C},1);
    Top_mean(1:40) = AveragedMatrix_excl(C,Top_Feats{C});

    Top_40{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];

end

% Indices of classifiers
wake = 0; N1 = 1; N2 = 3; N3 = 3; rem = 5;
CLASSIFIER = {[wake,N1] [wake,N2] [wake,N3] [wake,rem] [N1,N2] [N1,N3] [N1,rem] [N2,N3] [N2,rem] [N3,rem]};

% Which classifier do we want to plot?
ClNum = 10;  
Classif = Top_40{1,ClNum};  % Get the top features of the chosen classifier
Top_Feat = Top_Feats{1,ClNum};

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

% Get pairwise similarity matrix
numTopFeatures = 40;  % Number of top features
op_ind = Classif(:,3)'; % 3 = top_ID colum: indices of top 40 features (among the 7k)
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
clusterThreshold = 0.2;      % threshold at which split into clusters

% Ylabels
Top_mean = Classif(:,4)';  % Change it for later in the Ylabels section
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
           
                    
%% Get the top features for each classifier and each dataset
                   
% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
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


Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:12
    
    sub = Subs{D};
   
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))      

    for C = 1:10

        % Features reordering: from best to worst feature
        [~,I] = sort(Per_correct_mean_D_excl{1,D}(C,:)','descend');
        Top_Feats{C} = I(1:40);

        % Get the name and keyword associated with these features
        Top_name(1:40) = CodeString(Top_Feats{C},1);
        Top_key(1:40) = Keywords(Top_Feats{C},1);
        Top_ID(1:40) = YLabel(Top_Feats{C},1);
        Top_mean(1:40) = Per_correct_mean_D_excl{1,D}(C,Top_Feats{C});

        Top_40{D}{C} = [Top_name' Top_key' num2cell(Top_ID') num2cell(Top_mean')];

    end   
   
end
    

%% How many top features of a spe classifier are present across all datasets?

% load Matrix with TopFeat for each classifier and dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Top_Feat_Data_Class')  % all specifically removed featured across datasets (combined ('unique'))
whichClassif = 2;


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

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

load('HCTSA.mat', 'Operations')   % Simply to get the list of 7749 features 

equi_Top_Feat = setdiff(1:7749,spec_and_common_feat); % get the equivalent ranking of top feat after special-value features removed
Operations = Operations(equi_Top_Feat,:);  % get Operations of well-behaved features only

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';


% Convert 7749-based index (UniqueElem) to 5308-based index (Idx_UniqueElem)
for i = 1:length(UniqueElem)
    Idx_UniqueElem(i,:) = find(cell2mat(YLabel) == UniqueElem(i));
end

Top_Keywords = Keywords(Idx_UniqueElem(:,1));
Top_name = CodeString(Idx_UniqueElem(:,1));
Top_Keywords = [Top_name Top_Keywords num2cell(Top_Elem)];   % This top features may include features that are top in some datasets and special-value features in other datasets (that's why mean accuracy can be super low (around 7% if 11 datasets produced special values))


% Get % accuracy for each top feat (average of accuracy from ALL datasets)
Accuracy = [];
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')

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
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

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
    Top_Feats{C} = I(1:40);   
                     
    % Get the name and keyword associated with these features
    Top_name_aver(1:40) = CodeString(Top_Feats{C},1);
    Top_key_aver(1:40) = Keywords(Top_Feats{C},1);
    Top_ID_aver(1:40) = YLabel(Top_Feats{C},1);
    Top_mean_aver(1:40) = AveragedMatrix_excl(C,Top_Feats{C});

    Top_40_aver{C} = [Top_name_aver' Top_key_aver' Top_ID_aver' num2cell(Top_mean_aver')];

end


% Top_features of corresponding classifier; used in AverageMatrix
Top_aver = Top_40_aver{1,whichClassif};

% Find how many top features used in AverageMatrix are 
selected_feat_aver = find(ismember(cell2mat(Top_ID(1:40)),cell2mat(Top_aver(:,3))) == 1);
% ID_selected = Top_aver(selected_feat_aver,3);

% Plot distribution of top features across datasets
Featmode = cell2mat(Top_Keywords(:,4));   % y axis bar(Featmode)
Features = Top_Keywords(:,1);
Top_ID = Top_Keywords(:,3); 

h = bar(Featmode);

ax = gca;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:10:length(Features),'uni',0);
xticks(1:10:length(Features))
ylim([0 max(Featmode)+1])
xlabel('Top 40 features across all datasets (ranked from feature yielding highest to lowest accuracy) ')
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

g = bar(Fammode);

ax = gca;
ax.XTickLabels = arrayfun(@(a)string(a),Families(:,1:20),'uni',0);
xticks(1:length(Families))
xtickangle(35)
ylim([0 max(Fammode)+1])
xlabel('Top 20 families across all datasets (ranked from family yielding highest to lowest accuracy) ')
ylabel('Number of datasets sharing one top feature of the family')
title('Distribution of top families across datasets')

