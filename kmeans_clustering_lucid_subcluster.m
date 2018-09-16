lucid_base_configuration;

for SBJ_ID = config.subject_ids
    for MAIN_CLUSTERES_NUM=2:6

        for TARGET_CLUSTER_TO_SPLIT=1:MAIN_CLUSTERES_NUM    
            kmeans_clustering_configuration;

            BASE_PATH=strcat(config.base_dir, SBJ_ID, config.subject_secondary_id, config.run_base_folder);
            EDF_FILE=strcat(config.base_dir, 'edf_init_mat_files/', SBJ_ID, '_data.mat');
            STAGE_LOAD_FILENAME=strcat(BASE_PATH, filesep, 'HCTSA_N_', SBJ_ID, config.subject_secondary_id, ...
                '_1_EEG_Main_', num2str(MAIN_CLUSTERES_NUM), '_Clusters.mat');

            [feat_id, feat] = load_hctsa_reduced_ops(config.hctsa_reduced_ops_file, STAGE_LOAD_FILENAME);

            %% Use feat_id to select data from full op
            datamat = load(STAGE_LOAD_FILENAME,'TS_DataMat');
            orig_datamat = datamat.TS_DataMat;

            t = load(STAGE_LOAD_FILENAME, 'TimeSeries');
            ts = struct2table(t.TimeSeries);

            datamat = orig_datamat(C4_COL, feat_id);

            no_of_channels=config.no_of_channels_used;
            if no_of_channels == 1
                datamat = datamat;
            elseif no_of_channels == 2
                datamat = [datamat orig_datamat(EOG_COL,feat_id)];
            elseif no_of_channels == 3
                datamat = [datamat orig_datamat(EOG_COL,feat_id) orig_datamat(EMG_COL,feat_id)];        
            end

            idx=str2double(cellstr(table2array(ts(:, 2))));

            %% Lucid dreaming substage
            current_main_cluster = TARGET_CLUSTER_TO_SPLIT;
            current_idx = find(idx==current_main_cluster);
            sub_datamat=datamat(current_idx,:);

            [s, s_optk, ch, ch_optk, db, db_optk] = calculate_opt_k(sub_datamat);
            fprintf("====================================\n%s\n------------\n(Number of epochs=%d) \nSilhouette: %.04f (%d) \nCalinskiHarabasz: %.04f (%d) \nDaviesBouldin: %.04f (%d)\n\n", ...
                sprintf("Sub-cluster of main cluster %d", TARGET_CLUSTER_TO_SPLIT), ...
                size(sub_datamat, 1), ...
                s, ...
                s_optk, ...
                ch, ...
                ch_optk, ...
                db, ...
                db_optk);

            for sc = 1:size(config.sub_clusters_range, 2)
                subcluster=config.sub_clusters_range(sc);
                SUBSTAGE_SAVE_FILENAME=strcat('HCTSA_N_', SBJ_ID, config.subject_secondary_id, ...
                    '_TotalMain_', num2str(MAIN_CLUSTERES_NUM), ...
                    '_Cluster_', num2str(TARGET_CLUSTER_TO_SPLIT), '_1_EEG_', num2str(subcluster), ...
                    '_substages.mat');

                n_clust = subcluster;
                [idx, c, sse] = kmeans(sub_datamat,n_clust,'Distance','sqeuclidean',...
                                    'Display','off','Replicates',50,'MaxIter',500);    

                %% Save HTCSA_N trimmed (substages)
                t=load(hctsafile);
                dm = t.TS_DataMat;
                quality = t.TS_Quality;
                ts = t.TimeSeries;
                original_ts = struct2table(t.TimeSeries);

                single_channel_size = size(dm,1)/config.total_channels_in_hctsa_datamat;
                startIndex=1; % No trim

                % 1 EEG channel
                t.TS_DataMat = [dm(current_idx, feat_id)];
                t.TS_Quality = [quality(current_idx, feat_id)];
                t.TimeSeries = [ts(current_idx)];

                ts = t.TimeSeries;
                td=struct2table(ts);
                td(:,2)=cellstr([string(idx)]);
                ts=table2struct(td);

                t.TimeSeries = ts;

                ops=t.Operations;
                reduced_operations = [];

                feattable=struct2table(feat);
                for i = 1:size(feattable(:,1), 1)
                    f_name = table2array(feattable(i,2));
                    for j = 1:size(ops,1)
                        op = ops(j);
                        if string(op.Name) == string(cell2mat(f_name))
                           reduced_operations = [reduced_operations; op]; 
                        end
                    end
                end

                if (config.no_of_channels_used == 1)
                    t.Operations = [reduced_operations];
                elseif (config.no_of_channels_used == 3)
                    t.Operations = [reduced_operations reduced_operations reduced_operations];
                end

                save(strcat(BASE_PATH, filesep, SUBSTAGE_SAVE_FILENAME), '-struct', 't');
            end
        end
    end
    
    clearvars -except SBJ_ID config;
end