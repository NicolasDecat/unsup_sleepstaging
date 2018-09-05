for NUM_OF_MAIN_CLUSTERS=2:6

% homedir=pwd;
% %% Start HCTSA tools
% cd(HCTSA_DIR)
% startup
% cd(homedir)

% Comment the following two loops for single main clusters
for TARGET_MAIN_CLUSTER=1:NUM_OF_MAIN_CLUSTERS
for NUM_OF_SUB_CLUSTERS=2:5

%% Default configuration
kmeans_clustering_configuration;

%% Configuration

SBJ_ID='ME_N3';
SBJ_SECONDARY_ID='_bipolar';
FIGURE_VISIBLE = 'off';

HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
TARGET_FOLDER=strcat('/Volumes/Spaceship/Voss_Lucid/', SBJ_ID, SBJ_SECONDARY_ID, '/ALL_EEG');
NUM_SECONDS=30;
n_channels = 1;
% NUM_MAIN_CLUSTERS=5;
% NUM_OF_SUB_CLUSTERS=2;
% TARGET_MAIN_CLUSTER=4;

% MODE='MAIN_CLUSTER';
MODE='SUB_CLUSTER';
%MODE='PLOT_MAIN_CLUSTER';
%MODE='PLOT_SUB_CLUSTER';

if strcmp(MODE, 'MAIN_CLUSTER')
    fig_filename = sprintf('tsne_Main_%d_Clusters', NUM_OF_MAIN_CLUSTERS);
else
    fig_filename = sprintf('tsne_TotalMain_%d_Cluster_%d_SubCluster_%d', ...
        NUM_OF_MAIN_CLUSTERS, TARGET_MAIN_CLUSTER, NUM_OF_SUB_CLUSTERS);
end

TARGET_FILE=strcat('HCTSA_N_', SBJ_ID, SBJ_SECONDARY_ID, '_1_EEG_Main_', num2str(NUM_OF_MAIN_CLUSTERS), '_Clusters.mat');
if (strcmp(MODE, "SUB_CLUSTER"))
    TARGET_FILE=strcat('HCTSA_N_', SBJ_ID, SBJ_SECONDARY_ID, ...
        '_TotalMain_', num2str(NUM_OF_MAIN_CLUSTERS), ...
        '_Cluster_',num2str(TARGET_MAIN_CLUSTER), '_1_EEG_', num2str(NUM_OF_SUB_CLUSTERS), ...
        '_substages');
end
% STAGE_LOAD_FILENAME=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, '_Cluster_',num2str(TARGET_MAIN_CLUSTER), '_1_EEG_', num2str(NUM_OF_SUB_CLUSTERS), '_REM_substages.mat');

EDF_FILE=strcat('/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/', SBJ_ID, '_data.mat');

%% Main execution body
% distanceMetricRow = 'euclidean'; %
% linkageMethodRow = 'average'; %
% distanceMetricCol = 'corr_fast';
% linkageMethodCol = 'average'; %


%% Wanted operation names
OPS_FILE='reduced_ops.txt';
all_op = load(strcat(TARGET_FOLDER, filesep, TARGET_FILE),'Operations');
fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);

feat_name = features{1,2};

%% Check operation name, get feat_id
nn=0;
for n = 1:length(feat_name)
    op_name = char(feat_name(n));
    match_input_index = 0;
    match_test_index = 0;
    
    for i = 1:length(all_op.Operations)
        name = all_op.Operations(i).Name;
        if strcmp(op_name,name)
            match_input_index = i;
            break;
        end
    end
    
    if (match_input_index > 0)
        nn=nn+1;
        feat_id(nn) = match_input_index;
        feat(nn).id = match_input_index; % all_op.Operations(i).ID % Actual operation ID
        feat(nn).name = op_name;
    end
end
clear i n nn op_name name

%% Use feat_id to select data from full op
data = load(strcat(TARGET_FOLDER, filesep, TARGET_FILE));
orig_datamat = data.TS_DataMat;
ts = struct2table(data.TimeSeries);

if (strcmp(MODE, "SUB_CLUSTER"))
    datamat = orig_datamat;
else
    datamat = orig_datamat(C4_COL, feat_id);
end

if n_channels == 1
    datamat = datamat;
elseif n_channels == 2
    datamat = [datamat orig_datamat(EOG_COL,feat_id)];
elseif n_channels == 3
    datamat = [datamat orig_datamat(EOG_COL,feat_id) orig_datamat(EMG_COL,feat_id)];        
end

ts_labels = table2array(ts(:,1))';
ts_algo_groups = str2num(char(table2array(ts(:,2))));
unique_groups=unique(ts_algo_groups);  

clear cluster_indices;
for i = 1:length(unique_groups)
    cluster_indices{i} = find(ts_algo_groups==i);
end

max_elements = max(cellfun(@(x) size(x, 1), cluster_indices));

%% Loop through each clusters to gather standard deviation of each features
num_of_clusters = size(cluster_indices, 2);
% clusters_summary = struct([]);
% feature_table = struct2table(feat);
% 
% for g = 1:num_of_clusters
%     clusters_summary{g}.feature_std = std(datamat(cluster_indices{g}, :));
% end

%% Perform anova test among clusters
% features_diff = struct([]);
% for f = 1:size(feat, 2)
%     init_mat = NaN(max_elements, num_of_clusters);
%     
%     for g = 1:num_of_clusters
%         init_mat(1:size(cluster_indices{g}), g) = datamat(cluster_indices{g}, f);
%     end
%     
%     [p, tbl, stats] = anova1(init_mat, [1:num_of_clusters], 'off');
%     [c, m, h, nms] = multcompare(stats, 'Display', 'off');
%     
%     features_diff{f}.anova_p = p;
%     features_diff{f}.name = feat(f).name;
%     features_diff{f}.anova_F = tbl{2, 5};
%     features_diff{f}.pairwise_result = c;
% end

%% Analysis
% TOP_SIMILAR_FEATURE_WITHIN_CLUSTER = 10;
% ALPHA = 0.05;
% 
% for g = 1:num_of_clusters
%    fprintf("--------------- CLUSTER %d ---------------\n", g); 
%    
%    feature_std = clusters_summary{g}.feature_std;
%    [sort_values, sort_index] = sort(feature_std);
%    
%    for i = 1:TOP_SIMILAR_FEATURE_WITHIN_CLUSTER
%       % Check if the feature is unique (differ significantly among cluster)
%       % We are only interested in features that have significant diff
%       % in all pairwise tests.
%       feature_number = sort_index(i);
%       if features_diff{i}.anova_p < ALPHA
%           % There is at least one group differs significantly
%           % We checked if all pairwise are significant
%           if sum(features_diff{i}.pairwise_result(:, end) < 0.05) == ...
%                   size(features_diff{i}.pairwise_result, 1)
%              fprintf("%d: %s (std: %0.4f)\n", feature_number, ...
%                  features_diff{feature_number}.name,  sort_values(i));
%           end
%       end
%    end
% end

%% Plot two-dimensional analogues to 2D

Y = tsne(datamat, 'Algorithm', 'exact', 'Standardize', true, 'NumPCAComponents', 50);

set(0,'DefaultFigureVisible', FIGURE_VISIBLE);

%%
figure;
colors=num2cell(jet(length(unique_groups)* 5), 2);
cluster_colors = arrayfun(@(x) colors{x * 4}, [1:num_of_clusters], 'UniformOutput', false);

p=gscatter(Y(:,1),Y(:,2), ts_algo_groups, cell2mat(cluster_colors'), [], 20);

if (strcmp(MODE, "SUB_CLUSTER"))
    title(sprintf('Participant: %s - %d sub-clusters', SBJ_ID, NUM_OF_SUB_CLUSTERS), 'Interpreter', 'None');
else
    title(sprintf('Participant: %s - %d main clusters', SBJ_ID, NUM_OF_MAIN_CLUSTERS), 'Interpreter', 'None');
end

legend_labels = arrayfun(@(x) sprintf('Cluster %d', x), [1:num_of_clusters], 'UniformOutput', false);
legend(legend_labels);
set(gca, 'FontSize', 16);
grid on;
%% Plot grouping of data using the top features

saveas(gcf, strcat(TARGET_FOLDER, filesep, 'images', filesep, fig_filename), 'png');

% figure;
% %gscatter(datamat(:,165), datamat(:, 1), ts_algo_groups);
% feature_x = [120, 72, 165, 24, 59];
% feature_y = [77, 160, 113, 86, 54];
% feature_table=struct2table(feat);
% x_labels = num2str(strrep(table2array(feature_table(feature_x, 1)), '_', '\_'));
% y_labels = num2str(strrep(table2array(feature_table(feature_y, 1)), '_', '\_'));
% 
% set(gcf, 'Position', [200, 200, 1000, 800]);
% [h, ax, bigax] = gplotmatrix(datamat(:,feature_x), datamat(:, feature_y), ts_algo_groups, ...
%     'brgcm', '.x...', [], 'on', 'variable', x_labels, y_labels);
% set(ax(:, :), 'FontSize', 16);
% set(h, 'MarkerSize', 8);
% colormap(jet);

%% Plot 3-d scatter ploat

% figure;
% colors = brewermap(length(cluster_indices), 'Set1');
% 
% grid on;
% hold on;
% for k = 1:length(cluster_indices)
%     indices = cluster_indices{k};
%     %plot3(datamat(indices, 120), datamat
% end
% x = datamat(cluster_indices{1}, 120);
end
end
    
clear all;
end
