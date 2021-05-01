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

Subs = {'005'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
    Per_correct_mean = Per_correct_mean_D_excl{1,y(D)};

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
    subPerFig = 10; % subplots per figure
    numFeaturesDistr = 10;  % How many top features I want to plot

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
           if numel(Top_40{opi,1}) > 25
               Top_40{opi,1} = Top_40{opi,1}(1:25);
           end
            title({sprintf('[%u] %s: %1.1f%%',featHere(opi),string(Top_40{opi,1}),string(Top_40{opi,4}));...
                            ['(',Top_40{opi,2},')']},'interpreter','none')
        end
        
    end
    
     %%% Add line plot among subplots
    subplot(3,4,11)
    plot(Top_mean(1:10),'LineWidth',2)
    xlabel('Top 10 Features')
    ylabel('Accuracy (%)')

    ax.XTick = 1:10;
    ax.XTickLabels = arrayfun(@(a)num2str(a),0:10,'uni',0);
    labels = string(ax.XTickLabels); 
    ax.XTickLabels = labels; 
    ax.FontSize = 12; 
    title('Accuracy averaged across all classifiers')
    grid on
    
    %%% Alone
    figure; ax=gca;
    plot(Top_mean(1:10),'LineWidth',2)
    ylabel('Classificartion accuracy (%)')

    ax.XTick = 1:10;
    % xticklabels(string(Top_key(1:10)));
    ax.XTickLabel = {'wavelet','correlation','model','spectral','correlation','information','model','correlation','wavelet','forecasting'};
    xtickangle(45)
    ax.FontSize = 10; 
    grid on
    set(gca,'fontsize',20)
%     % Save
%     set(gcf,'color','white')
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
%     export_fig([fpath filesep 'Feat10'],'-r 300')

    
end


%%%%%%%% Plot the correlation matrix of top 10 features

numTopFeatures = numFeaturesDistr;  % Number of top features
ifeat = ifeat(1:10); % indices of top  features (among the 7k)

% Convert /7749 to /6606
load('HCTSA_N.mat','Operations')
for i = 1:numel(ifeat)
    op_ind(i) = find(Operations.ID == ifeat(i));
end

% Compute correlations based on hctsa responses
EEGonly = 1:size(TS_DataMat,1)/7;
TS_DataMat = TS_DataMat(EEGonly,op_ind);      % Only WB features and EEG channels 
Dij = BF_pdist(TS_DataMat','abscorr');  

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
                        'objectLabels',Top_name);
                    
set(gca,'fontsize',20)

    % Save
%     set(gcf,'color','white')
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
%     export_fig([fpath filesep 'Feat10_corrmat'],'-r 300')

%% Paper figure: violin plots of representative features

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

Subs = {'005'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
    % Per_correct_mean = iteration_svm_testing_accuracy_MEAN;  % if you want to plot top features from supervised clustering
    
    Per_correct_mean = Per_correct_mean_D_excl{1,y(D)};

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

    ifeat = Top_ID';

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot the correlation matrix of top 10 features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    numTopFeatures = 10;  % Number of top features
    ifeat = ifeat(1:10); % indices of top  features (among the 7k)

    % Convert /7749 to /6606
    load('HCTSA_N.mat','Operations')
    for i = 1:numel(ifeat)
        op_ind(i) = find(Operations.ID == ifeat(i));
    end

    % Compute correlations based on hctsa responses
    EEGonly = 1:size(TS_DataMat,1)/7;
    TS_DataMat = TS_DataMat(EEGonly,op_ind);      % Only WB features and EEG channels 
    Dij = BF_pdist(TS_DataMat','abscorr');  

    distanceMetric = 'abscorr';
    clusterThreshold = 0.2; % threshold at which split into clusters

    % Ylabels
    Top_mean = num2cell(mean(Top_40_Acc'));  % Change it for later in the Ylabels section
    YLabel = [];

    for i = 1:length(Top_ID)
        YLabel = [YLabel {sprintf('[%s] %s (%1.1f%%)',Top_key{i},Top_name{i},Top_mean{i})}];
    end

    % Plot
    [~,cluster_Groupi,RepFeat,ax1,ax2,f] = BF_ClusterDown_edited(Dij,'clusterThreshold',clusterThreshold,...
                            'whatDistance',distanceMetric,...
                            'objectLabels',YLabel);
      
   %%% Get position figure                     
   f.Position = [50 1 1300 900];   % [x y width height]
   ax1.Position = [0.77 0.1726 0.1023 0.7524];  % tree
   ax2.Position = [0.3 0.1726 0.5059 0.7524]; % matrix

   % save
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
%     export_fig([fpath filesep 'corrmat_10(439)'],'-r 300')
% 
%                         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot the violin plots of representative features
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ifeat = ifeat(RepFeat);  % ID of representative features (/7749)

    numClasses = 5; 
    subPerFig = length(RepFeat); % subplots per figure
    numFeaturesDistr = length(RepFeat);  % How many top features I want to plot

       
    % Set the colors to be assigned to groups:
    colors = GiveMeColors(numClasses);

    % Space the figures out properly:
    load('HCTSA_N.mat','TS_DataMat','TimeSeries','Operations')
    numFigs = ceil(numFeaturesDistr/subPerFig);

    % Make data structure for TS_SingleFeature
    data = struct('TS_DataMat',TS_DataMat,'TimeSeries',TimeSeries,...
                'Operations',Operations);


    % Get the indices of features to plot
    featHere = ifeat; 

    % Make the figure
    f = figure('color','w');
    f.Position(3:4) = [1353, 857];
    
    % Loop through features
    for opi = 1:length(featHere)
        subplot(ceil(length(featHere)/length(RepFeat)),length(RepFeat),opi);
        TS_SingleFeature_1D(data,featHere(opi),true,false);


%         title({sprintf('%s \n (%s) %1.1f%%',Top_name{RepFeat(opi)},Top_key{RepFeat(opi)},Top_mean{RepFeat(opi)});...
%                             },'interpreter','none')
       set(gca,'fontsize',20)
     
       title({sprintf('%s',Top_name{RepFeat(opi)})},'interpreter','none','fontsize',15)
                        
                               
    end


end

f.Position = [1 300 1430 400];   % [x y width height]

%Save
set(gcf,'color','white')
fpath = '/Users/nico/Documents/HCTSA/Analysis/violin';
export_fig([fpath filesep 'Feat10_violin'],'-r 300')


 %% Violin plots for all datasets
% 
% %%%%%% Step 1: from all 12 datasets, group the DataMat (feature values
% %%%%%% generated by top 40 features) by stage
% 
% %%% Get the ID of top 40 features for the 3 classifiers (N1, N2, REM)
% 
% % Compute top 40 features
% load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
% load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
% load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')
% 
% load HCTSA.mat   
% Operations_ID = [Operations.ID].';
% 
% equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); 
% Operations = Operations(equi_Top_Feat,:);                      
% 
% CodeString = {Operations.CodeString}.';  
% Keywords = {Operations.Keywords}.';
% YLabel = {Operations.ID}.';
% 
% TOP = 40;
% 
% for C = 1:10
%     
%     % Reorder features: from yielding theg highest to lowest accuracy
%     [~,I] = sort(AveragedMatrix_excl(C,:)','descend');
%     Top_Feats{C} = I(1:TOP);   
%                     
%     % Get the name and keyword associated with these features
%     Top_name(1:TOP) = CodeString(Top_Feats{C},1);
%     Top_key(1:TOP) = Keywords(Top_Feats{C},1);
%     Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
%     Top_mean(1:TOP) = AveragedMatrix_excl(C,Top_Feats{C});
% 
%     Top_40{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];
% end
% 
% 
% %%% Get the feature values (DataMat) of top 40 features in all 12 datasets
% 
% Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
% SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
% [~,y] = ismember(Subs,SUB);
% 
% for D = 1:length(Subs)  
%     
%     % ID of top 40 features for the 3 classifiers
%     N1N2 = cell2mat(Top_40{1,5}(:,3));   % ID out of 7749
%     N1REM = cell2mat(Top_40{1,7}(:,3));
%     N2REM = cell2mat(Top_40{1,9}(:,3));
% 
%  
%     sub = Subs{D};
%     
%     % Go to corresponding folder
%     cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
% 
%     % Load the normalized data and labels
%     load('HCTSA_N.mat','TimeSeries','TS_DataMat','Operations')
%     load(sprintf('ccshs_1800%s_annot.mat',sub), 'sleepstage')
%     length_EEG = 1:length(sleepstage);
%     
%     TS_DataMat = TS_DataMat(length_EEG,:);
%     TimeSeries = TimeSeries(length_EEG,:);
% 
%     % Convert feature ID: 7749 to number of WB features in a dataset (for
%     % loop to preserve order of ID)
%     for i = 1:40
%         N1N2(i,1) = find(ismember(Operations.ID,N1N2(i)));
%         N1REM(i,1) = find(ismember(Operations.ID,N1REM(i)));
%         N2REM(i,1) = find(ismember(Operations.ID,N2REM(i)));
%     end
%     
%     Classifiers = [N1N2 N1REM N2REM];
%     
%     % For each stage, get Datamat of top 40 features
%     numClasses = 5; 
%     classLabels = [0 1 2 3 5];
%     dataCell = cell(numClasses,1);
%     
%     for C = 1:3   % for the 3 classifiers
%         
%         for i = 1:numClasses
%             
%             % Row = stage, column = classifier, box = DataMat of top 40
%             % features for corresponding classifier
%             dataCell{i,C} = (TS_DataMat(sleepstage==classLabels(i),Classifiers(:,C)));
%         end
%         
%     end   
%       
%     
%     dataCell_all{D} = dataCell;
%     clear N1N2 N1REM N2REM dataCell
% end
% 
% 
% % Group epochs from all datasets, per stage
% all_wake = []; all_N1 = []; all_N2 = []; all_N3 = []; all_rem = [];
% 
% for C = 1:3
%     
%     wake = []; N1 = []; N2 = []; N3 = []; rem = [];
% 
%      for D = 1:12
%     
%         wake = [wake;dataCell_all{1,D}(1,C)];
%         N1 = [N1;dataCell_all{1,D}(2,C)];
%         N2 = [N2;dataCell_all{1,D}(3,C)];
%         N3 = [N3;dataCell_all{1,D}(4,C)];
%         rem = [rem;dataCell_all{1,D}(5,C)];
% 
%      end
%      
%      % Epochs from all datasets gathered
%      all_wake{C} = wake;
%      all_N1{C} = N1;
%      all_N2{C} = N2;
%      all_N3{C} = N3;
%      all_rem{C} = rem;
% 
% end
% 
% % Re-arrange the epochs and average
% wake = []; N1 = []; N2 = []; N3 = []; rem = [];
% 
% for C = 1:3
%     
%     for D = 1:12
%         
%     wake = [wake;all_wake{1,C}{D,1}];
%     N1 = [N1;all_N1{1,C}{D,1}];
%     N2 = [N2;all_N2{1,C}{D,1}];
%     N3 = [N3;all_N3{1,C}{D,1}];
%     rem = [rem;all_rem{1,C}{D,1}];
%     
%     end
%     
%     % For each classifer, you have the feature values from all datasets 
%     % grouped by stage
%     ALL_wake{C} = wake;
%     ALL_N1{C} = N1;
%     ALL_N2{C} = N2;
%     ALL_N3{C} = N3;
%     ALL_rem{C} = rem;
%     
% end
% 
% % 3 columns, each is a classifer. 5 rows for each stage
% dataCell = [];
% dataCell = [ALL_wake;ALL_N1;ALL_N2;ALL_N3;ALL_rem];
% 
% 
% %%%%%%%% Step 2: plot violin plots
% 
% Classif = 1;  % N1N2
% Classif_ID = [5 7 9];
% 
% numClasses = 5; 
% subPerFig = 10; % subplots per figure
% numFeaturesDistr = 10;  % How many top features I want to plot
% 
% ifeat = Classifiers(:,Classif);
% 
% % Set the colors to be assigned to groups:
% colors = GiveMeColors(numClasses);
% 
% % Space the figures out properly:
% numFigs = ceil(numFeaturesDistr/subPerFig);
% 
% for figi = 1:numFigs
%     
%     if figi*subPerFig > TOP
%         break % We've exceeded number of features
%     end
%     
%     % Get the indices of features to plot
%     r = ((figi-1)*subPerFig+1:figi*subPerFig);
%     
%     if figi==numFigs % filter down for last one
%         r = r(r<=numFeaturesDistr);
%     end
%     
%     featHere = ifeat(r); % features to plot on this figure
%     feature = figi;
%     
%     f = figure('color','w');
%     f.Position(3:4) = [1353, 857];
%     
%     % Loop through features
%     for opi = 1:length(featHere)
%         subplot(ceil(length(featHere)/4),4,opi);
%         TS_SingleFeature_allD(dataCell,featHere(opi),true,false,feature,C);
% %        if numel(Top_40{opi,1}) > 25
% %            Top_40{opi,1} = Top_40{opi,1}(1:25);
% %        end
% %         title({sprintf('[%u] %s: %1.1f%%',featHere(opi),string(Top_40{opi,1}),string(Top_40{opi,4}));...
% %                         ['(',Top_40{opi,2},')']},'interpreter','none')
%     end
% 
% end
% 
% 
% 
% 

%% Line plot for top 10 for each classifier 

% Go to corresponding current folder
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
Dataset = 3;
sub = SUB{Dataset};
    
cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
%%%%%%  Prepare top 10
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')
load HCTSA.mat   
Operations_ID = [Operations.ID].';

equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); % remove SV features
Operations = Operations(equi_Top_Feat,:);                      % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

TOP = 40;

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')

CL = {'Wake vs N1', 'Wake vs N2', 'Wake vs N3', 'Wake vs REM', 'N1 vs N2', 'N1 vs N3', 'N1 vs REM', 'N2 vs N3', 'N2 vs REM', 'N3 vs REM'};

figure;    
[ha, pos] = tight_subplot(2,5,[.115 .05],[.1 .05],[.05 .05]);

TOP = 10;

%%%%% Select a classifier 
for C = 1:10
    
    Per_correct_mean = Per_correct_mean_D_excl{1,Dataset}(C,:);

    [~,I] = sort(Per_correct_mean','descend');
    Top_Feats{C} = I(1:TOP);   
                    
    % Get the name and keyword associated with these features
    Top_name(1:TOP) = CodeString(Top_Feats{C},1);
    Top_key(1:TOP) = Keywords(Top_Feats{C},1);
    Top_ID(1:TOP) = YLabel(Top_Feats{C},1);
    Top_mean(1:TOP) = Per_correct_mean(Top_Feats{C});

    Top_10{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];

    %%%%% Line plot of top 10 features
    axes(ha(C)); 
    
    plot(Top_mean,'LineWidth',2)
    xlabel('Top 10 Features')
    ylabel('Accuracy (%)')

    ax.XTick = 0:10;
    ax.XTickLabels = arrayfun(@(a)num2str(a),0:10,'uni',0);
    labels = string(ax.XTickLabels); 
    ax.XTickLabels = labels; 
    ax.FontSize = 12; 
    title(CL(C))
    grid on

    
     %%% Violin plot

    numClasses = 5; 
    subPerFig = 10; % subplots per figure
    numFeaturesDistr = 10;  % How many top features I want to plot

    ifeat = Top_ID';

    % Set the colors to be assigned to groups:
    colors = GiveMeColors(numClasses);

    % Space the figures out properly:
    numFigs = ceil(numFeaturesDistr/subPerFig);

    load HCTSA_N
    
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
        featHere = cell2mat(ifeat(r)); % features to plot on this figure
        % featHere = find(ismember(YLabel, featHere));  % Make equivalent
        % Make the figure
        f = figure('color','w');
        f.Position(3:4) = [1353, 857];
        % Loop through features
        for opi = 1:length(featHere)
            NameFeat = Top_10{1,C}{opi,1};
            subplot(ceil(length(featHere)/4),4,opi);
            TS_SingleFeature_1D(data,featHere(opi),true,false);
           
           if numel(NameFeat) > 25
               NameFeat = NameFeat(1:25);
           end
            title({sprintf('[%u] %s: %1.1f%%',featHere(opi),string(NameFeat),string(Top_10{1,C}{opi,4}));...
                            ['(',Top_10{1,C}{opi,2},')']},'interpreter','none')
        end
        
    end
    
    %%% Add line plot among subplots
    subplot(3,4,11)
    plot(Top_mean,'LineWidth',2)
    xlabel('Top 10 Features')
    ylabel('Accuracy (%)')

    ax.XTick = 0:10;
    ax.XTickLabels = arrayfun(@(a)num2str(a),0:10,'uni',0);
    labels = string(ax.XTickLabels); 
    ax.XTickLabels = labels; 
    ax.FontSize = 12; 
    title(CL(C))
    grid on
    
    % Save
    sgtitle(sprintf('Top 10 Features for: %s',CL{C}))
    fpath = '/Users/nico/Documents/HCTSA/Analysis/violin/Dataset439';
%     saveas(gca,fullfile(fpath,sprintf('top10_%s(%s)',string(CL{C}),sub)),'jpg')
    
    %%%%%%%% Add correlation matrix among subplots

    numTopFeatures = numFeaturesDistr;  % Number of top features
    ifeat = cell2mat(ifeat(1:10)); % indices of top  features (among the 7k)

    % Convert /7749 to /6606
    load('HCTSA_N.mat','Operations')
    for i = 1:numel(ifeat)
        op_ind(i) = find(Operations.ID == ifeat(i));
    end

    % Compute correlations based on hctsa responses
    EEGonly = 1:size(TS_DataMat,1)/7;
    TS_DataMat = TS_DataMat(EEGonly,op_ind);      % Only WB features and EEG channels 
    Dij = BF_pdist(TS_DataMat','abscorr');  

    distanceMetric = 'abscorr';
    clusterThreshold = 0.2; % threshold at which split into clusters

    % Ylabels
    Top_Mean = Top_10{1,C}(:,4)';  % Change it for later in the Ylabels section
    Ylabel = [];

    for i = 1:length(Top_ID)
        Ylabel = [Ylabel {sprintf('[%s] %s (%1.1f%%)',Top_key{i},Top_name{i},Top_Mean{i})}];
    end

    % Plot
    [~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                            'whatDistance',distanceMetric,...
                            'objectLabels',Ylabel);
    title(sprintf('Dependencies between %u top features (%u clusters), %s',...
                            numTopFeatures,length(cluster_Groupi),CL{C}))
         
                                      
end




