

%% Turn plotting on
set(0,'DefaultFigureVisible','on')

%% Standard Plotting: accuracy matrix (try things here)

figure;
imagesc(AveragedMatrix_excl)
% title(sprintf('Classification performance per feature (Dataset %s)',sub));
% title('Classification performance per feature');

ax = gca;
ax.XTick = 1:500:6056;
ax.YTick = 1:10;
set(ax, 'TickLength', [0 0]);
% ax.XTickLabels = strseq('f',1:100)';
ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:6006,'uni',0);
ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
ylabel('Binary classifiers');
xlabel('features');
ax.XAxisLocation = 'bottom';
ax.FontSize = 12; 

colormap 'default'
colorbar

% Save figure
exportgraphics(gca,'/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/AveragedMatrix_HQ.jpg','Resolution',200)


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

% Get the name of all 7749 features 
all_op = load('HCTSA.mat','Operations');
all_op_name = {all_op.Operations.CodeString}.';
    
for D = 1:length(Subs)   
    
    sub = Subs{D};

    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))  
    
    % Get the name of all retained operations (exclude SV features)
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

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/common_features_removed')
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

% Original features - names
all_op = load('HCTSA.mat','Operations');
all_op_name = {all_op.Operations.CodeString}.';
all_7749 = 1:7749;

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Configuration
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))  
    
    retained_op = load('HCTSA_N.mat','Operations');
    retained_op_name = retained_op.Operations.CodeString;

    all_op_idx = retained_op.Operations{:,4};                             
    removed_feat_idx{D} =  all_7749(~ismember(all_7749,all_op_idx));   
    
    specifically_removed{D} = setdiff(removed_feat_idx{1,D},val);  % Get only features that are not commonly found across datasets (the "specific ones")

end

%% Get all features removed (common + spe across datasets)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/common_features_removed')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/specifically_removed_features')

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
spec_and_common_feat = sort([val All_spe]);


%% Reconstruct the matrix by including the removed special-value features (for one dataset)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/allfeat_removed')  % All features removed
Dataset = 12;
spec_and_common_feat = removed_feat_idx{1,Dataset};

InsertCol = zeros(10,1);  % column that will be inserted

for x = 1:length(spec_and_common_feat)
    
    Part1 = Per_correct_mean(:,1:spec_and_common_feat(x)-1);
    Part2 = [InsertCol Per_correct_mean(:,spec_and_common_feat(x):end)];  % Add the 'zeros' column (act as a special-value feature)
    Per_correct_mean = ([Part1 Part2]);  % then merge the parts, with the zeros column in sandwich between the 2
    
end

%% Include in the accuracy matrix the specifically removed features (blue bars)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/allfeat_removed')  % All features removed
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/common_features_removed')   % Only features commonly removed in all datasets
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/specifically_removed_features')  % Only features specifically removed in eachd dataset

InsertCol = zeros(10,1);

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))

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

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146).mat') 
InsertCol = zeros(10,1);

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))

    for x = 1:length(removed_feat_idx{1,D})

        Part1 = Per_correct_mean(:,1:removed_feat_idx{1,D}(x)-1);
        Part2 = [InsertCol Per_correct_mean(:,removed_feat_idx{1,D}(x):end)];
        Per_correct_mean = ([Part1 Part2]);
        removed_feat_idx{1,D} = removed_feat_idx{1,D};

    end
    
    % Remove both the specific and commonly removed features
    Per_correct_mean(:,spec_and_common_feat) = [];   
    Per_correct_mean_D{D} = Per_correct_mean;
    
end


% Stack all 12 matrices along 3rd dimension
StackedMatrix = cat(3,Per_correct_mean_D{1,1},Per_correct_mean_D{1,2},Per_correct_mean_D{1,3},Per_correct_mean_D{1,4},Per_correct_mean_D{1,5},Per_correct_mean_D{1,6},Per_correct_mean_D{1,7},Per_correct_mean_D{1,8},Per_correct_mean_D{1,9},Per_correct_mean_D{1,10},Per_correct_mean_D{1,11},Per_correct_mean_D{1,12});
AveragedMatrix = mean(StackedMatrix,3);

%%% 5603 features are kept in total (7779 minus the 2146 special value
%%% features)

% Then plot AveragedMatrix (w/o special value features; preserved ranking)
figure;
imagesc(AveragedMatrix_excl)
title('Mean classification performance across datasets');

ax = gca;
ax.XTick = 1:500:5603;
ax.YTick = 1:10;
% ax.XTickLabels = strseq('f',1:100)';
ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:5603,'uni',0);
ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
ylabel('Binary classifiers');
xlabel('features');
ax.XAxisLocation = 'bottom';

colormap 'default'
colorbar


%% Reconstruct each dataset to a 7749-matrix and plot all features (special values in red)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/All_unique_specificity_feat_combined')  % all specifically removed featured across datasets (combined ('unique'))
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/allfeat_removed') % For each of the 12 cells, All features removed for the corresponding dataset) 
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/common_features_removed')   % For each of the 12 cells, only features commonly removed (shared by all datasets) for the corresponding dataset
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/specifically_removed_features')  % For each of the 12 cells, only features specifically removed in the corresponding dataset

InsertCol = zeros(10,1);

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)  
    
    sub = Subs{D};
    
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean(Dataset %s)',sub))

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
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')  % all specifically removed featured across datasets (combined ('unique'))
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


%% Data Matrices in 5603 (no SV features)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl') 

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:12
    
    sub = Subs{D};

    Per_correct_mean = Per_correct_mean_D_excl{1,D};
    
    figure;
    imagesc(Per_correct_mean)
    title(sprintf('Classification performance per feature (Dataset %s)',sub));
    ax = gca;
    ax.XTick = 1:500:5603;
    ax.YTick = 1:10;
    ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:5603,'uni',0);
    ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
    ylabel('Binary classifiers');
    xlabel('features');
    ax.XAxisLocation = 'bottom';
    colormap 'default'
    colorbar

    % Save figure
    fpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Figure_accuracy_per_feat/AccMatrix_noSV_5603';
    saveas(gca,fullfile(fpath,sprintf('AccperFeat(100)_5603_%s_crossval',sub)),'fig')
    saveas(gca,fullfile(fpath,sprintf('AccperFeat(100)_5603_%s_crossval',sub)),'jpg')
 
    
end


%% Plots unsupervised / supervised / per-feature / all-feature accuracies

%% Line plots %%%%%%%%%%%%%%%%%%

% Accuracies dataset 001
uns_all = flip([73.0 77.1 92.1 74.5 63.2 93.8 59.9 84.5 63.8 91.1])';
uns_one = flip([54.2 56.3 64.0 56.1 58.3 69.4 54.5 62.5 54.2 69.0])';
sup_all = flip([89.1 92.3 98.2 95.5 86.4 99.5 82.5 93.2 84.1 99.1])';
sup_one = flip([58.4 59.9 69.7 59.8 57.5 69.7 54.4 64.0 55.6 66.4])';

% Sort from best to worst classifier (based on unsup all)
[~,I] = sort(uns_all,'ascend');
uns_all = uns_all(I);
uns_one = uns_one(I);
sup_all = sup_all(I);
sup_one = sup_one(I);

% Line plot
figure; 
h = plot(1:10,uns_all,'LineWidth',1.3,'Color',[0.3010 0.7450 0.9330]);  % light blue
hold on
i = plot(1:10,sup_all,'LineWidth',1.3,'Color',[0 0.4470 0.7410]);  % dark blue
hold on
j = plot(1:10, uns_one,'LineWidth',1.3,'Color',[0.9350 0.580 0.3840]);  % light red
hold on
k = plot(1:10,sup_one,'LineWidth',1.3,'Color',[0.6350 0.0780 0.1840]);  % dark red
hold off

% legend('Unsup - using all features','SVM - using all features','Unsup - one feature at a time','SVM - one feature at a time','Location','eastoutside')
xlabel('Classifiers')
ylabel('Classification accuracy (%)')

ax = gca;
ax.XTick = 1:11;
ax.XTickLabels = {'N3 vs REM','N2 vs REM','N3 vs N2','N1 vs REM','N3 vs N1','N1 vs N2','W vs REM','N3 vs W','W vs N2','W vs N1',''};
ax.XTickLabels = ax.XTickLabels(I);
xtickangle(30)
ylim([45 100])
yline(50,'--','chance level');
ax.FontSize = 14;

%% Rainbow plots %%%%%%%%%%%%
 
% Load data
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean(Dataset 001)')

% Acc per feature, Wake vs N3 classifier  
Data{1,1} = Per_correct_mean(3,:); 
% Acc per feature, Wake vs N1 classifier 
Data{2,1} = Per_correct_mean(1,:); 
% Acc per feature, REM vs N3 classifier  
Data{3,1} = Per_correct_mean(10,:); 
% Acc per feature, REM vs N2 classifier 
Data{4,1} = Per_correct_mean(9,:); 


%%%% Plot Rainbow graph

% Set colors
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);
cl(3, :) = cb(2, :);

figure; 

%%%%%% Left plot: W vs N1

subplot(2,2,1);  

% Wake vs N1: accuracy per feature
l1 = raincloud_plot(Data{2,1}, 'box_on', 1, 'color', cb(7,:), 'alpha', 0.5,...    
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);

yl = ylim;  % used later to get default max on y axis

hold on
m = bar(70.7,yl(2));   % Bar for all-features accuracy (unsup)
hold on 
n = bar(89.1,yl(2));   % Bar for all-features accuracy (sup)
hold on
mean2 = mean(Data{2,1});
b = xline(mean2,'Color',cb(7,:),'LineWidth',2.5);

title('Wake vs N1')
o = line([50 50],[yl(1) yl(2)],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2);  % Line for 50% accuracy
ylim([-0.015 yl(2)]);
legend([l1{1} b m(1) n(1) o(1)], {'Accuracy per feature (unsup)','Mean accuracy per feature (unsup)','Accuracy all features (unsup)','Accuracy all features (SVM)','50% accuracy'},'Location','northwest','FontSize',10);
xlabel('Accuracy')
ylabel('Distribution')

%%%%% Right plot: Wake vs N3

subplot(2,2,2);  

% Wake vs N3: accuracy per feature
h1 = raincloud_plot(Data{1,1}, 'box_on', 1, 'color', cb(7,:), 'alpha', 0.5,...    
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);

yl = ylim;  % used later to get default max on y axis

hold on
i = bar(90.6,yl(2));   % Bar for all-features accuracy (unsup)
hold on 
j = bar(98.2,yl(2));   % Bar for all-features accuracy (sup)
hold on
mean1 = mean(Data{1,1});
a = xline(mean1,'Color',cb(7,:),'LineWidth',2.5);
ylim([0 1])

title('Wake vs N3')
k = line([50 50],[yl(1) yl(2)],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2);  % Line for 50% accuracy
ylim([-0.009 yl(2)]);
xlabel('Accuracy')
ylabel('Distribution')


%%%%%% Bottom Left plot: N2 vs REM

subplot(2,2,3);  

% Wake vs N1: accuracy per feature
ll1 = raincloud_plot(Data{4,1}, 'box_on', 1, 'color', cb(7,:), 'alpha', 0.5,...    
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);

yl = ylim;  % used later to get default max on y axis

hold on
mm = bar(62.3,yl(2));   % Bar for all-features accuracy (unsup)
hold on 
nn = bar(84.1,yl(2));   % Bar for all-features accuracy (sup)
hold on
mean4 = mean(Data{4,1});
c = xline(mean4,'Color',cb(7,:),'LineWidth',2.5);

title('REM vs N2')
oo = line([50 50],[yl(1) yl(2)],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2);  % Line for 50% accuracy
ylim([-0.02 yl(2)]);
xlabel('Accuracy')
ylabel('Distribution')

%%%%% Bottom Right plot: REM vs N3

subplot(2,2,4);   

% REM vs N3: accuracy per feature
hh1 = raincloud_plot(Data{3,1}, 'box_on', 1, 'color', cb(7,:), 'alpha', 0.5,...    
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);

yl = ylim;  % used later to get default max on y axis

hold on
ii = bar(93.1,yl(2));   % Bar for all-features accuracy (unsup)
hold on 
jj = bar(99.1,yl(2));   % Bar for all-features accuracy (sup)
hold on
mean3 = mean(Data{3,1});
d = xline(mean3,'Color',cb(7,:),'LineWidth',2.5);

title('REM vs N3')
kk = line([50 50],[yl(1) yl(2)],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2);  % Line for 50% accuracy
ylim([-0.008 yl(2)]);
xlabel('Accuracy')
ylabel('Distribution')


%% Fix 870

TimeSeries = struct2table(TimeSeries);
for i = 1:length(TimeSeries.Data)
    TimeSeries.Data(i,1) = TimeSeries.Data{i,1};
end
TimeSeries.Data = cell2mat(TimeSeries.Data);

%% Plot data matrix of accuracy across datasets for given classifier

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl') 
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

NumFeat = 5603;
whichClassif = 5;

for D = 1:12
    
    sub = Subs{D};
    Per_correct_mean(D,:) = Per_correct_mean_D_excl{1,D}(whichClassif,1:NumFeat);  % e.g., N1 vs N2, 100 features
    
end

figure;
imagesc(Per_correct_mean)
title('Classification performance per feature across datasets');
ax = gca;
ax.XTick = 1:50:NumFeat;
ax.YTick = 1:12;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:50:NumFeat,'uni',0);
ax.YTickLabels = {'001','005','439','458','596','748','749','752','604','807','821','870'};
ylabel('Datasets');
xlabel('features');
ax.XAxisLocation = 'bottom';
colormap 'default'
colorbar

%% Reorder the figure above: most consistent features across datasets first

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl') 
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

NumFeat = 5603;
whichClassif = 5;

for D = 1:12
    
    sub = Subs{D};
    Per_correct_mean(D,:) = Per_correct_mean_D_excl{1,D}(whichClassif,1:NumFeat);  % e.g., N1 vs N2, 100 features

end

% get mean and standard deviation across datasets for given classifer
MEAN = mean(Per_correct_mean);
STD = std(Per_correct_mean);

% Compute consistency and sort from most to least consistent
Consistent = MEAN./STD;
[~,I] = sort((Consistent)','descend');
Consistent = Consistent(:,I);
Per_correct_mean = Per_correct_mean(:,I);


% Plot
figure;
imagesc(Consistent)
title('Classification performance per feature across datasets');
ax = gca;
ax.XTick = 1:50:NumFeat;
ax.YTick = 1:12;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:50:NumFeat,'uni',0);
ax.YTickLabels = {'001','005','439','458','596','748','749','752','604','807','821','870'};
ylabel('Datasets');
xlabel('features');
ax.XAxisLocation = 'bottom';
colormap 'default'
colorbar
