lucid_base_configuration;

for SBJ_ID = config.subject_ids
    for MODE = ["MAIN", "SUB_CLUSTER"]
    for NUM_OF_MAIN_CLUSTERS=2:6

        [MAX_TARGET_MAIN_CLUSTER, MAX_SUB_CLUSTERS] = determine_max_clusters(MODE, NUM_OF_MAIN_CLUSTERS);

        % Comment the following two loops for single main clusters
        for TARGET_MAIN_CLUSTER=1:MAX_TARGET_MAIN_CLUSTER
            for NUM_OF_SUB_CLUSTERS=2:MAX_SUB_CLUSTERS
        
            kmeans_clustering_configuration;

            %% Configuration
            TARGET_FOLDER=strcat(config.base_dir, SBJ_ID, config.subject_secondary_id, config.run_base_folder);

            if strcmp(MODE, 'MAIN')
                fig_filename = sprintf('tsne_Main_%d_Clusters', NUM_OF_MAIN_CLUSTERS);
            else
                fig_filename = sprintf('tsne_TotalMain_%d_Cluster_%d_SubCluster_%d', ...
                    NUM_OF_MAIN_CLUSTERS, TARGET_MAIN_CLUSTER, NUM_OF_SUB_CLUSTERS);
            end

            TARGET_FILE=strcat('HCTSA_N_', SBJ_ID, config.subject_secondary_id, '_1_EEG_Main_', num2str(NUM_OF_MAIN_CLUSTERS), '_Clusters.mat');
            if (strcmp(MODE, "SUB_CLUSTER"))
                TARGET_FILE=strcat('HCTSA_N_', SBJ_ID, config.subject_secondary_id, ...
                    '_TotalMain_', num2str(NUM_OF_MAIN_CLUSTERS), ...
                    '_Cluster_',num2str(TARGET_MAIN_CLUSTER), '_1_EEG_', num2str(NUM_OF_SUB_CLUSTERS), ...
                    '_substages');
            end
            
            [feat_id, feat] = load_hctsa_reduced_ops(config.hctsa_reduced_ops_file, ...
                strcat(TARGET_FOLDER, filesep, TARGET_FILE));

            %% Use feat_id to select data from full op
            data = load(strcat(TARGET_FOLDER, filesep, TARGET_FILE));
            orig_datamat = data.TS_DataMat;
            ts = struct2table(data.TimeSeries);

            if size(ts, 1) > size(orig_datamat, 1)
               ts = ts(1:size(orig_datamat, 1), :);
            end

            if config.no_of_channels_used == 1
                feat_id = [feat_id];
            elseif config.no_of_channels_used == 3
                feat_id = [feat_id, length(feat_id)+feat_id, (length(feat_id)*2)+feat_id];
            end

            if (strcmp(MODE, "SUB_CLUSTER"))
                datamat = orig_datamat;
            else
                datamat = orig_datamat(C4_COL, feat_id);
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

            set(0,'DefaultFigureVisible', 'off');              

            %%
            figure;
            %colors=num2cell(jet(length(unique_groups)* 5), 2);
            colors=GiveMeColors(length(unique_groups));

            cluster_colors = arrayfun(@(x) colors{x}, [1:num_of_clusters], 'UniformOutput', false);

            p=gscatter(Y(:,1),Y(:,2), ts_algo_groups, cell2mat(cluster_colors'), [], 20);

            if (strcmp(MODE, "SUB_CLUSTER"))
                title(sprintf('Participant: %s - %d sub-clusters', SBJ_ID, NUM_OF_SUB_CLUSTERS), 'Interpreter', 'None');
            else
                title(sprintf('Participant: %s - %d main clusters', SBJ_ID, NUM_OF_MAIN_CLUSTERS), 'Interpreter', 'None');
            end

            legend_labels = arrayfun(@(x) sprintf('Cluster %d', x), [1:num_of_clusters], 'UniformOutput', false);
            legend(legend_labels);
            set(gca, 'FontSize', 16);

            set(gca,'xtick',[])
            set(gca,'ytick',[])
            %grid on;

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
        
    clearvars -except SBJ_ID MODE config;
    end
    end
end

function [MAX_TARGET_MAIN_CLUSTER, MAX_SUB_CLUSTERS] = determine_max_clusters(mode, NUM_OF_MAIN_CLUSTERS)
    if strcmp(mode, "MAIN")
        MAX_TARGET_MAIN_CLUSTER=1;
        MAX_SUB_CLUSTERS = 2;
    else
        MAX_TARGET_MAIN_CLUSTER=NUM_OF_MAIN_CLUSTERS;
        MAX_SUB_CLUSTERS = 7;
    end
end
