%% Violin Plots: distribution of top features across classes

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Top_Feat_Data_Class')  

% Load TS_DataMat 
load('HCTSA_N.mat','Operations','TS_DataMat','TimeSeries')  

% Get equivalent indices (only WB features)
equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features

% Trim Operations to WB features
EEGonly = 1:size(TS_DataMat,1);
Operations = Operations(Idx_WB_Feat,:);

% Average TS_DataMat across all 12 datasets
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

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

TS_DataMat = TS_DataMat(1:1088,Idx_WB_Feat);
    
CodeString = Operations.CodeString;  
Keywords = Operations.Keywords;
YLabel = Operations.ID;

% Features reordering: from best to worst feature
means = mean(AveragedMatrix_excl);  
[~,I] = sort((means)','descend');
AveragedMatrix = AveragedMatrix_excl(:,I);

% Get best 40 features 
Top_Feat = I(1:40); 

Top_name(1:40) = CodeString(Top_Feat,1);
Top_key(1:40) = Keywords(Top_Feat,1);
Top_ID(1:40) = YLabel(Top_Feat,1);

Top_40 = [Top_name' Top_key' num2cell(Top_ID')];


% Parameters plot
numClasses = 5; 
subPerFig = 16; % subplots per figure
numFeaturesDistr = 16;  % How many top features I want to plot

ifeat = Top_ID';

% Set the colors to be assigned to groups:
colors = GiveMeColors(numClasses);

% Space the figures out properly:
numFigs = ceil(numFeaturesDistr/subPerFig);

% Make data structure for TS_SingleFeature
data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
            'Operations',Operations);

for figi = 1:numFigs
    if figi*subPerFig > length(ifeat)
        break % We've exceeded number of features
    end
    % Get the indices of features to plot
    r = ((figi-1)*subPerFig+1:figi*subPerFig);
    if figi==numFigs % filter down for last one
        r = r(r<=numFeaturesDistr);
    end
    featHere = ifeat(r); % features to plot on this figure
    featHere = find(ismember(cell2mat(YLabel), featHere));  % Make equivalent
    % Make the figure
    f = figure('color','w');
    f.Position(3:4) = [1353, 857];
    % Loop through features
    for opi = 1:length(featHere)
        subplot(ceil(length(featHere)/4),4,opi);
        TS_SingleFeature(data,Operations.ID(featHere(opi),:),true,false);
    end
end


%% Violin plots for one dataset

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'439'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

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


%%%  Parameters plot
%TS_DataMat = TS_DataMat(1:1374,:);

numClasses = 5; 
subPerFig = 16; % subplots per figure
numFeaturesDistr = 16;  % How many top features I want to plot

ifeat = Top_ID';

% Set the colors to be assigned to groups:
colors = GiveMeColors(numClasses);

% Space the figures out properly:
numFigs = ceil(numFeaturesDistr/subPerFig);

% Make data structure for TS_SingleFeature
data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
            'Operations',Operations);

for figi = 1:numFigs
    if figi*subPerFig > length(ifeat)
        break % We've exceeded number of features
    end
    % Get the indices of features to plot
    r = ((figi-1)*subPerFig+1:figi*subPerFig);
    if figi==numFigs % filter down for last one
        r = r(r<=numFeaturesDistr);
    end
    featHere = ifeat(r); % features to plot on this figure
    % featHere = find(ismember(YLabel, featHere));  % Make equivalent
    % Make the figure
    f = figure('color','w');
    f.Position(3:4) = [1353, 857];
    % Loop through features
    for opi = 1:length(featHere)
        subplot(ceil(length(featHere)/4),4,opi);
        TS_SingleFeature_1D(data,featHere(opi),true,false);
        title({sprintf('[%u] %s: %s',featHere(opi),string(Top_40{opi,1}),string(Top_40{opi,4}));...
                        ['(',Top_40{opi,2},')']},'interpreter','none')

    end
end

%% Violin plots for one dataset, one classifier

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Top_Feat_Data_Class')
load HCTSA_N.mat

Dataset = 1;
Classif = 10;
Top40 = Top_40{1,Dataset}{1,Classif};

%%%  Parameters plot
%TS_DataMat = TS_DataMat(1:1374,:);

numClasses = 5; 
subPerFig = 16; % subplots per figure
numFeaturesDistr = 16;  % How many top features I want to plot

ifeat = Top40(:,3)';

% Set the colors to be assigned to groups:
colors = GiveMeColors(numClasses);

% Space the figures out properly:
numFigs = ceil(numFeaturesDistr/subPerFig);

% Make data structure for TS_SingleFeature
data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
            'Operations',Operations);

for figi = 1:numFigs
    if figi*subPerFig > length(ifeat)
        break % We've exceeded number of features
    end
    % Get the indices of features to plot
    r = ((figi-1)*subPerFig+1:figi*subPerFig);
    if figi==numFigs % filter down for last one
        r = r(r<=numFeaturesDistr);
    end
    featHere = ifeat(r); % features to plot on this figure
    % featHere = find(ismember(YLabel, featHere));  % Make equivalent
    % Make the figure
    f = figure('color','w');
    f.Position(3:4) = [1353, 857];
    % Loop through features
    for opi = 1:length(featHere)
        subplot(ceil(length(featHere)/4),4,opi);
        TS_SingleFeature(data,featHere(opi),true,false);
        title({sprintf('[%u] %s: %s',Top40{opi},Top40{opi,1},Top40{opi,4});...
                        ['(',Top40{opi,2},')']},'interpreter','none')

    end
end