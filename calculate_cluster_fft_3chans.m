lucid_base_configuration;

for SBJ_ID = config.subject_ids
    for MODE = ["MAIN", "SUB_CLUSTER"]

        EDF_FILE=strcat(config.base_dir, 'edf_init_mat_files/', SBJ_ID, '_data.mat');
        load(EDF_FILE);

        for NUM_OF_MAIN_CLUSTERS=2:6

            [MAX_TARGET_MAIN_CLUSTER, MAX_SUB_CLUSTERS] = determine_max_clusters(MODE, NUM_OF_MAIN_CLUSTERS);

            for TARGET_MAIN_CLUSTER=1:MAX_TARGET_MAIN_CLUSTER
                for NUM_OF_SUB_CLUSTERS=2:MAX_SUB_CLUSTERS
    
                kmeans_clustering_configuration;

                if strcmp(MODE, "MAIN")
                    fig_filename = sprintf('Main_%d_Clusters', NUM_OF_MAIN_CLUSTERS);
                else
                    fig_filename = sprintf('TotalMain_%d_Cluster_%d_SubCluster_%d', ...
                        NUM_OF_MAIN_CLUSTERS, TARGET_MAIN_CLUSTER, NUM_OF_SUB_CLUSTERS);
                end
                
                BASE_PATH=strcat(config.base_dir, SBJ_ID, config.subject_secondary_id, config.run_base_folder);

                if (strcmp(MODE, "MAIN"))
                    MAIN_CLUSTER_FILE_PREFIX=strcat(BASE_PATH, filesep, 'HCTSA_N_', SBJ_ID, config.subject_secondary_id, '_1_EEG_Main_', ...
                        num2str(NUM_OF_MAIN_CLUSTERS), '_Clusters');
                    STAGE_LOAD_FILENAME=strcat(MAIN_CLUSTER_FILE_PREFIX,  '.mat');
                    STAGE_LOAD_CSV=strcat(MAIN_CLUSTER_FILE_PREFIX,  '.csv');
                else
                    SUB_CLUSTER_FILE_PREFIX=strcat(BASE_PATH, filesep, 'HCTSA_N_', SBJ_ID, config.subject_secondary_id, ...
                        '_TotalMain_', num2str(NUM_OF_MAIN_CLUSTERS), ...
                        '_Cluster_',num2str(TARGET_MAIN_CLUSTER), '_1_EEG_', num2str(NUM_OF_SUB_CLUSTERS), ...
                        '_substages');
                    STAGE_LOAD_FILENAME=strcat(SUB_CLUSTER_FILE_PREFIX, '.mat');
                    STAGE_LOAD_CSV=strcat(SUB_CLUSTER_FILE_PREFIX, '.csv');
                end

                kmeans_clustering_configuration_channels;
                plot_power_spectrum(STAGE_LOAD_FILENAME, header, data, epoch_seconds, ...
                    timeseries_sampling_rate, chans, STAGE_LOAD_CSV, 0);
                saveas(gcf, strcat(BASE_PATH, filesep, 'images', filesep, fig_filename), 'png');

                end
            end
        end
        
        clearvars -except SBJ_ID config MODE data header;
        
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