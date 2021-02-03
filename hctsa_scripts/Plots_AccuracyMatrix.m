

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


%% Plots unsupervised / supervised / per-feature / all-feature accuracies

%%% Line plots




