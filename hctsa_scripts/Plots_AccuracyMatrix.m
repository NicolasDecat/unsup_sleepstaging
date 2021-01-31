

%% Turn plotting on
set(0,'DefaultFigureVisible','on')

%% Standard Plotting: accuracy matrix (try things here)

figure;
imagesc(Per_correct_mean)
% title(sprintf('Classification performance per feature (Dataset %s)',sub));
title('Classification performance per feature (Dataset 870)');

ax = gca;
ax.XTick = 1:500:5892;
ax.YTick = 1:10;
% ax.XTickLabels = strseq('f',1:100)';
ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:5892,'uni',0);
ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
ylabel('Binary classifiers');
xlabel('features');
ax.XAxisLocation = 'bottom';

colormap 'default'
colorbar


%% Line plot of accuracy for each binary classifier across features

%%%%% All binary classifiers in the same plot

%%% Average every 100 features to make the plot more readable

for classifier = 1:10
    x = Per_correct_mean(classifier,:)';
    S = numel(x);
    xx = reshape(x(1:S - mod(S, 100)), 100, []);
    y(:,classifier)  = sum(xx, 1).' / 100;
end
y = y';  

%%% Then plot y

figure; plot(y(:,:)','LineWidth',1.3)
% figure; plot(AUC_per_feature(:,:)')
xlim([0 60])
legend('W vs N1','W vs N2','W vs N3','W vs REM','N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM','Location','eastoutside')
xticklabels({'0','1000','2000','3000','4000','5000','6000'})
xlabel('features')
ylabel('accuracy')


%%%%% One binary classifier only in the plot

figure; plot(Per_correct_mean(1,:))
xlim([0 6006])
ylim([0 100])
legend('Wake vs N1');   % W vs N1','W vs N2','W vs N3','W vs REM','N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM')
xlabel('features');
ylabel('accuracy');


%% Plot average accuracy across datasets (wrong, can't overlap those because feature ordering not the same acorss datasets) 

%%%  All 12 datasets must have the same length as the shortest dataset (752) 

Per_correct_mean001 = Per_correct_mean001(:,1:5858);
Per_correct_mean005 = Per_correct_mean005(:,1:5858);
Per_correct_mean439 = Per_correct_mean439(:,1:5858);
Per_correct_mean458 = Per_correct_mean458(:,1:5858);
Per_correct_mean596 = Per_correct_mean596(:,1:5858);
Per_correct_mean604 = Per_correct_mean604(:,1:5858);
Per_correct_mean748 = Per_correct_mean748(:,1:5858);
Per_correct_mean749 = Per_correct_mean749(:,1:5858);
Per_correct_mean752 = Per_correct_mean752(:,1:5858);
Per_correct_mean807 = Per_correct_mean807(:,1:5858);
Per_correct_mean821 = Per_correct_mean821(:,1:5858);
Per_correct_mean870 = Per_correct_mean870(:,1:5858);

%%% Stack all 12 matrices along a 3rd dimension

StackedMatrix = cat(3,Per_correct_mean001,Per_correct_mean005,Per_correct_mean439,Per_correct_mean458,Per_correct_mean596,Per_correct_mean604,Per_correct_mean748,Per_correct_mean749,Per_correct_mean752,Per_correct_mean807,Per_correct_mean821,Per_correct_mean870);
AveragedMatrix = mean(StackedMatrix,3);

%%% Then plot the averaged matrix

figure;
imagesc(AveragedMatrix)
title('Mean classification performance across datasets');

ax = gca;
ax.XTick = 1:500:5858;
ax.YTick = 1:10;
% ax.XTickLabels = strseq('f',1:100)';
ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:5858,'uni',0);
ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
ylabel('Binary classifiers');
xlabel('features');
ax.XAxisLocation = 'bottom';

colormap 'default'
colorbar

%% Features reordering based on TS_CLUSTER (op_clust)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D')

TS_Cluster
Per_correct_mean = Per_correct_mean_D{1,1};
Per_correct_mean = Per_correct_mean(:,(op_clust.ord)');


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
   

%% Reorder the features to get proper alignment 

%%%%% How: put back all special value features: now all feat ordered. Then
%%%%% delete all those special value features for all datasets

%% Identify the features removed by TS_Quality

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Get the name of all 7749 features 
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))  
    all_op = load('HCTSA.mat','Operations');
    all_op_name = {all_op.Operations.CodeString}.';
    
    % Get the name of all retained operations (exclude special value
    % features)
    retained_op = load('HCTSA_N.mat','Operations');
    retained_op_name = retained_op.Operations.CodeString;
    
    % Get the ID of these operations (just in case)
    ID_retained_feat = retained_op.Operations.ID;     % Useful to get feature names of TS_Cluster accuracy matrix

    all_op_idx = retained_op.Operations{:,4};   % Index of retained feat
    all_7749 = 1:7749;                          % Index of all features
    
    % For each dataset, get the index of special-value features removed 
    removed_feat_idx{D} =  all_7749(~ismember(all_7749,all_op_idx));   

end

% Get the features that were removed across all datasets ("common removed
% features")  = val
Z = removed_feat_idx;
[val]=intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(intersect(Z{1},Z{2}),Z{3}),Z{4}),Z{5}),Z{6}),Z{7}),Z{8}),Z{9}),Z{10}),Z{11}),Z{12});

%%% Mean number of features removed per dataset: 1775.8
%%% 75% of removed features are common to all datasets. About 400 features were removed in each dataset and not in others. 


%% Identify the features removed specifically in single dataset

%%% Can run this section only after getting removed_feat_idx (common 
%%% removed features  (run section above)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Configuration
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))  
    all_op = load('HCTSA.mat','Operations');
    all_op_name = {all_op.Operations.CodeString}.';
    
    retained_op = load('HCTSA_N.mat','Operations');
    retained_op_name = retained_op.Operations.CodeString;

    all_op_idx = retained_op.Operations{:,4};   
    all_7749 = 1:7749;                          
    removed_feat_idx{D} =  all_7749(~ismember(all_7749,all_op_idx));   
    
    specifically_removed{D} = setdiff(removed_feat_idx{1,D},val);  % Get only features that are not commonly found across datasets (the "specific ones")

end

%% Get all features removed (common + spe across datasets)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Configuration
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))
    
    ALLfeat_removed{1,D} = sort([val specifically_removed{1,D}]);
end


%%% list of all special value features removed
ALLSPE = [];

for D = 1:12
    ALLSPE = [ALLSPE specifically_removed{1,D}];
end

All_spe = unique(ALLSPE);

% List of all features removed across datasets
spec_and_common_feat2 = sort([val All_spe]);


%% Reconstruct the matrix by including the removed special-value features (for one dataset)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')

InsertCol = zeros(10,1);  % column that will be inserted

for x = 1:length(removed_feat_idx{:,:})
    
    Part1 = Per_correct_mean(:,1:removed_feat_idx{1,1}(x)-1);
    Part2 = [InsertCol Per_correct_mean(:,removed_feat_idx{1,1}(x):end)];  % Add the 'zeros' column (act as a special-value feature)
    Per_correct_mean = ([Part1 Part2]);  % then merge the parts, with the zeros column in sandwich between the 2
    
end


%% Include in the accuracy matrix the specifically removed features (blue bars)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/allfeat_removed')  % All features removed
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')   % Only features commonly removed in all datasets
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')  % Only features specifically removed in eachd dataset

InsertCol = zeros(10,1);

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))

    for x = 1:length(removed_feat_idx{1,D})

        Part1 = Per_correct_mean(:,1:removed_feat_idx{1,D}(x)-1);
        Part2 = [InsertCol Per_correct_mean(:,removed_feat_idx{1,D}(x):end)];
        Per_correct_mean = ([Part1 Part2]);
        removed_feat_idx{1,D} = removed_feat_idx{1,D};

    end
    
    Per_correct_mean(:,val) = [];   
    Per_correct_mean_D{D} = Per_correct_mean;  % 12 cells for each dataset: accuracy matrix displaying the specifically removed features (blue bars) (commonly removed features are not displayed)
    
end

%%% If want to display the commonly-removed features, don't do "Per_correct_mean(:,val) = [];" above   


%% Remove both common and specifically removed features (preserved ranking)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset

InsertCol = zeros(10,1);

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))

    for x = 1:length(removed_feat_idx{1,D})

        Part1 = Per_correct_mean(:,1:removed_feat_idx{1,D}(x)-1);
        Part2 = [InsertCol Per_correct_mean(:,removed_feat_idx{1,D}(x):end)];
        Per_correct_mean = ([Part1 Part2]);
        removed_feat_idx{1,D} = removed_feat_idx{1,D};

    end
    
    % Remove both the specific and commonly removed features
    spec_and_common_feat = [val All_spec_feat];
    Per_correct_mean(:,spec_and_common_feat) = [];   
    Per_correct_mean_D{D} = Per_correct_mean;
    
end


% Stack all 12 matrices along 3rd dimension
StackedMatrix = cat(3,Per_correct_mean_D{1,1},Per_correct_mean_D{1,2},Per_correct_mean_D{1,3},Per_correct_mean_D{1,4},Per_correct_mean_D{1,5},Per_correct_mean_D{1,6},Per_correct_mean_D{1,7},Per_correct_mean_D{1,8},Per_correct_mean_D{1,9},Per_correct_mean_D{1,10},Per_correct_mean_D{1,11},Per_correct_mean_D{1,12});
AveragedMatrix = mean(StackedMatrix,3);

%%% 5308 features are kept in total (7774 minus the 2441 special value
%%% features)

% Then plot AveragedMatrix (w/o special value features; preserved ranking)
figure;
imagesc(AveragedMatrix)
title('Mean classification performance across datasets');

ax = gca;
ax.XTick = 1:500:5308;
ax.YTick = 1:10;
% ax.XTickLabels = strseq('f',1:100)';
ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:5308,'uni',0);
ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
ylabel('Binary classifiers');
xlabel('features');
ax.XAxisLocation = 'bottom';

colormap 'default'
colorbar


%% Reconstruct each dataset to a 7749-matrix and plot all features (special values in red)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset

InsertCol = zeros(10,1);

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))

    % I'll include both the commonly and specifically removed features for
    % this dataset
    spec_common = sort([val specifically_removed{1,D}]);
    
    for x = 1:length(spec_common)

        Part1 = Per_correct_mean(:,1:spec_common(x)-1);
        Part2 = [InsertCol Per_correct_mean(:,spec_common(x):end)];
        Per_correct_mean = ([Part1 Part2]);

    end
    
    % Store
    Per_correct_mean_D{D} = Per_correct_mean;
    
end


% Stack all 12 matrices along 3rd dimension
StackedMatrix = cat(3,Per_correct_mean_D{1,1},Per_correct_mean_D{1,2},Per_correct_mean_D{1,3},Per_correct_mean_D{1,4},Per_correct_mean_D{1,5},Per_correct_mean_D{1,6},Per_correct_mean_D{1,7},Per_correct_mean_D{1,8},Per_correct_mean_D{1,9},Per_correct_mean_D{1,10},Per_correct_mean_D{1,11},Per_correct_mean_D{1,12});
AveragedMatrix = mean(StackedMatrix,3);


%%%% %All special value features will be marked as red
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')  % all specifically removed featured across datasets (combined ('unique'))
AveragedMatrix(:,spec_and_common_feat) = 0;   % (give value 0 to be marked in red by colormap) / give value [] to get matrix without special value features


% Then plot 
figure;
imagesc(AveragedMatrix)
title('Mean classification performance across datasets');

ax = gca;
ax.XTick = 1:500:7749;
ax.YTick = 1:10;
% ax.XTickLabels = strseq('f',1:100)';
ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:7749,'uni',0);
ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
ylabel('Binary classifiers');
xlabel('features');
ax.XAxisLocation = 'bottom';

cmap = colormap;
cmap(1,1:3) = [1 0 0];   % common special-value features in red
colormap(cmap)           % Activate it
colorbar



%% Return the top features from averaged matrix 

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_incl_all_feat_removed(7749)_all_datasets')

% Features reordering: from best to worst feature
means = mean(AveragedMatrix);  % AveragedMatrix = 7749-long including special value features
[~,I] = sort((means)','descend');
AveragedMatrix = AveragedMatrix(:,I);

% Get best X features 
Top_Feat = I(1:500);   % Get top 10 best features

% Get the name and keyword associated with these features
load('HCTSA.mat', 'Operations')   % Simply to get the list of 7749 features 

CodeString = {Operations.CodeString}.';   % Get name
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

Top_name = CodeString(Top_Feat,1);
Top_key = Keywords(Top_Feat,1);
Top_ID = YLabel(Top_Feat,1);
Top_mean = mean(AveragedMatrix(:,1:500))';

Top_10 = [Top_name Top_key Top_ID num2cell(Top_mean)];


%% Same as above but top features for each classifier

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_incl_all_feat_removed(7749)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D')

load('HCTSA.mat', 'Operations')   % Simply to get the list of 7749 features 
CodeString = {Operations.CodeString}.';   % Get name
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';


%%% Top 3 per classifier (over datasets)

for C = 1:10
  
% Features reordering: from best to worst feature
[~,I] = sort(AveragedMatrix(C,:)','descend');
Top_Feat = I(1:3);   % Get the best feature
        
% Get the name and keyword associated with these features
Top_name(1:3) = CodeString(Top_Feat,1);
Top_key(1:3) = Keywords(Top_Feat,1);
Top_ID(1:3) = YLabel(Top_Feat,1);
Top_mean(1:3) = AveragedMatrix(C,Top_Feat);

Top_3{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];

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

%% on getting Per_correct_mean_D

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
    Per_correct_mean_D{D} = Per_correct_mean;
    
end


%% Top 40 features for one Dataset (only WB Features)

% Matrices with only WB features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

Subs = {'604'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

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

%%%%% Get the AveragedMatrix with all 12 datasets, 7749 features, and their
%%%%% top 40 features

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_incl_all_feat_removed(7749)_all_datasets')

% Features reordering: from best to worst feature
means = mean(AveragedMatrix);  % AveragedMatrix = 7749-long including special value features
[~,I] = sort((means)','descend');
AveragedMatrix = AveragedMatrix(:,I);

% Get best 40 features 
Top_Feat = I(1:40);  

% Get the name and keyword associated with these features
load('HCTSA.mat', 'Operations')   % Simply to get the list of 7749 features 

% Get the name and keyword associated with these features
CodeString = {Operations.CodeString}.';   
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

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
    load('HCTSA.mat', 'TS_DataMat')    
    
    % Store
    DataMat{D} = TS_DataMat;
    
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
op_ind = Top_ID'; % indices of top 40 features (among the 7k)

op_ind = cell2mat(op_ind);

% Compute correlations based on hctsa responses
Dij = BF_pdist(AveragedDataMat(:,op_ind)','abscorr'); 
% EEGonly = size(AveragedDataMat,1)/7;
% Dij = BF_pdist(TS_DataMat(1:EEGonly,op_ind)','abscorr');   

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

%%%%% Plot the corr matrix between hctsa values of 2 stages of a classifier

% Copy pasted script to get 7749-matrix
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset


Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once
% (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'};
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
       
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D')
    
end

Per_correct_mean = Per_correct_mean_D{1,y};

% Load Data
cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
load HCTSA.mat

% Get top features per classifier

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

for C = 1:10

    % Features reordering: from best to worst feature
    [~,I] = sort(Per_correct_mean(C,:)','descend');
    Top_Feat = I(1:40);

    % Get the name and keyword associated with these features
    Top_name(1:40) = CodeString(Top_Feat,1);
    Top_key(1:40) = Keywords(Top_Feat,1);
    Top_ID(1:40) = YLabel(Top_Feat,1);
    Top_mean(1:40) = Per_correct_mean(C,Top_Feat);

    Top_40{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];

end

% Indices of classifiers
wake = 0; N1 = 1; N2 = 3; N3 = 3; rem = 5;
CLASSIFIER = {[wake,N1] [wake,N2] [wake,N3] [wake,rem] [N1,N2] [N1,N3] [N1,rem] [N2,N3] [N2,rem] [N3,rem]};

% Which classifier do we want to plot?
ClNum = 3;  
Classif = Top_40{1,ClNum};  % Get the top features of the chosen classifier

% hctsa values of this classif
load(sprintf('ccshs_1800%s_annot.mat',sub), 'sleepstage')   
stage1 = find(sleepstage == CLASSIFIER{ClNum}(1));  % indices of epochs that are labeled stage1
stage2 = find(sleepstage == CLASSIFIER{ClNum}(2));  % indices of epochs that are labeled stage2

DataMat = sort([stage1;stage2]);  % indices of epochs stage1 + stage2

%%%%%%%% Plot the correlation matrix

% Get pairwise similarity matrix
numTopFeatures = 40;  % Number of top features
op_ind = Classif(:,3)'; % 3 = top_ID colum: indices of top 40 features (among the 7k)

op_ind = cell2mat(op_ind);

% Compute correlations based on hctsa responses
DataMat = DataMat';
EEGonly = 1:length(TimeSeries)/7;
% Dij = BF_pdist(TS_DataMat(DataMat,op_ind)','abscorr');    % Only EEG channels that are epochs of the stages
Dij = BF_pdist(TS_DataMat(EEGonly,op_ind)','abscorr');    % all EEG channels

% Dij = BF_pdist(TS_DataMat(1:EEGonly,op_ind)','abscorr');   

distanceMetric = 'abscorr';
clusterThreshold = 0.2; % threshold at which split into clusters

% Ylabels
Top_mean = Classif(:,4)';  % Change it for later in the Ylabels section
YLabel = [];
for i = 1:length(Top_ID)
    YLabel = [YLabel {sprintf('[%s] %s (%1.1f%%)',Classif{i,2},Classif{i,1},Classif{i,4})}];
end

% Plot
[~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',YLabel);
title(sprintf('Dependencies between %u top features (organized into %u clusters)',...
                        numTopFeatures,length(cluster_Groupi)))
           
                    
%% Corr matrix averaged over all 12 datasets on specific classifiers

%%%%% Get the AveragedMatrix with all 12 datasets between 2 spe stages 

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_incl_all_feat_removed(7749)_all_datasets')

% Load Data
load HCTSA.mat

% Get top features per classifier

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

for C = 1:10

    % Features reordering: from best to worst feature
    [~,I] = sort(AveragedMatrix(C,:)','descend');
    Top_Feat = I(1:40);

    % Get the name and keyword associated with these features
    Top_name(1:40) = CodeString(Top_Feat,1);
    Top_key(1:40) = Keywords(Top_Feat,1);
    Top_ID(1:40) = YLabel(Top_Feat,1);
    Top_mean(1:40) = AveragedMatrix(C,Top_Feat);

    Top_40{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];

end

% Indices of classifiers
wake = 0; N1 = 1; N2 = 3; N3 = 3; rem = 5;
CLASSIFIER = {[wake,N1] [wake,N2] [wake,N3] [wake,rem] [N1,N2] [N1,N3] [N1,rem] [N2,N3] [N2,rem] [N3,rem]};

% Which classifier do we want to plot?
ClNum = 5;  % Wake vs N1
Classif = Top_40{1,ClNum};  % Get the top features of the chosen classifier


%%%% Now get the hctsa values averaged over all datasets

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
SUB = {'001','005','439','458','596','748','749','752','604','807','821','870'
    };
[~,y] = ismember(Subs,SUB);

for D = 1:length(Subs)  
 
    sub = Subs{D};
    
    % Go to corresponding current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Load TS_DataMat 
    load('HCTSA.mat', 'TS_DataMat')   % Simply to get the list of 7749 features 
    
    % Store
    DataMat{D} = TS_DataMat;
    
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

% Compute correlations based on hctsa responses (can't do otherwise than
% chosing all epochs, because epochs N1 N2 etc across datasets are diff.
EEGonly = 1:size(AveragedDataMat,1)/7;  
Dij = BF_pdist(AveragedDataMat(EEGonly,op_ind)','abscorr');    % Only EEG channels that are epochs of the stages

distanceMetric = 'abscorr';
clusterThreshold = 0.2; % threshold at which split into clusters

% Ylabels
Top_mean = Classif(:,4)';  % Change it for later in the Ylabels section
YLabel = [];
for i = 1:length(Top_ID)
    YLabel = [YLabel {sprintf('[%s] %s (%1.1f%%)',Classif{i,2},Classif{i,1},Classif{i,4})}];
end

% Plot
[~,cluster_Groupi] = BF_ClusterDown(Dij,'clusterThreshold',clusterThreshold,...
                        'whatDistance',distanceMetric,...
                        'objectLabels',YLabel);
title(sprintf('Dependencies between %u top features (organized into %u clusters)',...
                        numTopFeatures,length(cluster_Groupi)))
           
                    
%% Get the top features for each classifier and each dataset
                   
% Copy pasted script to get 7749-matrix
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset

InsertCol = zeros(10,1);

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% This is just to obtain the right index when not all Subs are ran at once (use y)
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
    Per_correct_mean_D{D} = Per_correct_mean;
    
end


Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:12
    
    sub = Subs{D};
   
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))      % Load Data
    load HCTSA.mat

    % Get top features per classifier
    CodeString = {Operations.CodeString}.';  
    Keywords = {Operations.Keywords}.';
    YLabel = {Operations.ID}.';

    for C = 1:10

        % Features reordering: from best to worst feature
        [~,I] = sort(Per_correct_mean_D{1,D}(C,:)','descend');
        Top_Feat = I(1:40);

        % Get the name and keyword associated with these features
        Top_name(1:40) = CodeString(Top_Feat,1);
        Top_key(1:40) = Keywords(Top_Feat,1);
        Top_ID(1:40) = YLabel(Top_Feat,1);
        Top_mean(1:40) = Per_correct_mean_D{1,D}(C,Top_Feat);

        Top_40{D}{C} = [Top_name' Top_key' Top_ID' num2cell(Top_mean')];

    end   
   
end
    

%% How many top features of a spe classifier are present across all datasets?

% load Matrix with TopFeat for each classifier and dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Top_Feat_Data_Class')  % all specifically removed featured across datasets (combined ('unique'))
whichClassif = 1;


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

% Get the keywords assigned to Feat ID
load('HCTSA.mat','Operations');
Keywords = {Operations.Keywords}.';
CodeString = {Operations.CodeString}.';

Top_Keywords = Keywords(Top_Elem(:,1));
Top_name = CodeString(Top_Elem(:,1));
Top_Keywords = [Top_name Top_Keywords num2cell(Top_Elem)];   % This top features may include features that are top in some datasets and special-value features in other datasets (that's why mean accuracy can be super low (around 7% if 11 datasets produced special values))


% Get % accuracy for each top feat (average of accuracy from ALL datasets)
Accuracy = [];
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D')

for ID = 1:length(UniqueElem)
    
    for D = 1:12
        Accuracy(D) = Per_correct_mean_D{1,D}(whichClassif,UniqueElem(ID));
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

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_incl_all_feat_removed(7749)_all_datasets')

% Load Data
load HCTSA.mat

% Get top features per classifier

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

for C = 1:10

    % Features reordering: from best to worst feature
    [~,I] = sort(AveragedMatrix(C,:)','descend');
    Top_Feat_aver = I(1:40);

    % Get the name and keyword associated with these features
    Top_name_aver(1:40) = CodeString(Top_Feat_aver,1);
    Top_key_aver(1:40) = Keywords(Top_Feat_aver,1);
    Top_ID_aver(1:40) = YLabel(Top_Feat_aver,1);
    Top_mean_aver(1:40) = AveragedMatrix(C,Top_Feat_aver);

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


%% Try on HCTSA_N TS_DataMat


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
    
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
  
    Per_correct_mean = Per_correct_mean_D_excl{D};
    
end

% load AveragedMatrix without special value features, and the list of
% special value features
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5308)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/ALL_removed_feat(2441)')

load('HCTSA_N.mat', 'Operations')    

% From Operations, take only well behaved features
equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % get the equivalent ranking of top feat after special-value features removed

Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of well-behaved features

Operations = Operations(Idx_WB_Feat,:);  % get Operations of well-behaved features only

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';


% Reorder features
means = mean(Per_correct_mean);
[~,I] = sort((means)','descend');
Per_correct_mean = Per_correct_mean(:,I);
    
% Get the best 40 features
Top_Feat = I(1:40);   

% Get the name and keyword associated with these features
CodeString = Operations.CodeString;   % Get name
Keywords = Operations.Keywords;
YLabel = Operations.ID;

Top_name(1:40) = CodeString(Top_Feat,1);
Top_key(1:40) = Keywords(Top_Feat,1);
Top_ID(1:40) = YLabel(Top_Feat,1);

Top_40 = [Top_name' Top_key' num2cell(Top_ID')];

% Get mean over classifiers
for F = 1:length(Top_40)
    for C = 1:10
        Top_40_Acc(F,C) = mean(Per_correct_mean(C,F));
    end
end

Top_mean = mean(Top_40_Acc');

Top_40 = [Top_40 num2cell(Top_mean')];


%%%%%%%% Plot the correlation matrix

%%% Run section above (Top 40 features with dataset 001)

% Get pairwise similarity matrix
numTopFeatures = 40;  % Number of top features
op_ind = Top_ID'; % indices of top 40 features (among the 7k)

% Compute correlations based on hctsa responses
load('HCTSA_N.mat', 'TS_DataMat')
EEGonly = 1:size(TS_DataMat,1)/7;    % Only EEG channels
TS_DataMat = TS_DataMat(EEGonly,Idx_WB_Feat);   % Only WB features 

Dij = BF_pdist(TS_DataMat(:,Top_Feat)','abscorr');    
% Dij = BF_pdist(TS_DataMat(1:EEGonly,op_ind)','abscorr');   

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
                    

  
 