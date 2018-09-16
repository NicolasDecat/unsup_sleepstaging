lucid_base_configuration;

for SBJ_ID = config.subject_ids
   
    opt_k_report_str = "";

    for n_clust = 2:6

        %% CONFIGURATION
        n_channels = config.no_of_channels_used;

        kmeans_clustering_configuration;

        STAGE_SAVE_FILENAME=strcat('HCTSA_N_', SBJ_ID, config.subject_secondary_id, '_1_EEG_Main_', num2str(n_clust), '_Clusters.mat');
        EDF_FILE=strcat(config.base_dir, 'edf_init_mat_files/', SBJ_ID, '_data.mat');
        BASE_OUTPUT_FOLDER=strcat(config.base_dir, SBJ_ID, config.subject_secondary_id, config.run_base_folder);

        [feat_id, feat] = load_hctsa_reduced_ops(config.hctsa_reduced_ops_file, hctsafile);

        datamat = load(hctsafile,'TS_DataMat');
        orig_datamat = datamat.TS_DataMat;

        datamat = orig_datamat(C4_COL, feat_id);

        if n_channels == 1
            datamat = datamat;
        elseif n_channels == 2
            datamat = [datamat orig_datamat(EOG_COL,feat_id)];
        elseif n_channels == 3
            datamat = [datamat orig_datamat(EOG_COL,feat_id) orig_datamat(EMG_COL,feat_id)];        
        end

        [s, s_optk, ch, ch_optk, db, db_optk] = calculate_opt_k(datamat);
        opt_k_report_str = strcat(opt_k_report_str, ...
            sprintf("====================================\n%s\n------------\n(Number of epochs=%d) \nSilhouette: %.04f (%d) \nCalinskiHarabasz: %.04f (%d) \nDaviesBouldin: %.04f (%d)\n\n", ...
                "Main cluster: ", ...
                size(datamat, 1), ...
                s, ...
                s_optk, ...
                ch, ...
                ch_optk, ...
                db, ...
                db_optk));

        %% Perform k-means clustering
        [idx, c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                            'Display','off','Replicates',50,'MaxIter',500);

        %% Save HTCSA_N trimmed for HCTSA analysis (assigned labels)
        t=load(hctsafile);
        dm = t.TS_DataMat;
        quality = t.TS_Quality;
        ts = t.TimeSeries;

        single_channel_size = size(dm,1)/no_of_channels;

        startIndex=1; % No trim
        endIndex=single_channel_size;

        if n_channels == 1
            %1 EEG channel
            t.TS_DataMat = [dm(startIndex:endIndex, feat_id)];
            t.TS_Quality = [quality(startIndex:endIndex, feat_id)];
            t.TimeSeries = [ts(startIndex:endIndex)];
        elseif n_channels == 3
            %3 channels
            t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(single_channel_size+startIndex:single_channel_size*2, feat_id) dm(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
            t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(single_channel_size+startIndex:single_channel_size*2, feat_id) quality(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
            t.TimeSeries = [ts(startIndex:single_channel_size);ts(single_channel_size+startIndex:single_channel_size*2);ts(single_channel_size*2+startIndex:single_channel_size*3)];
        end

        ts = t.TimeSeries;
        td=struct2table(ts);

        if n_channels == 1
            td(:,2)=cellstr([string(idx)]);
        elseif n_channels == 3
            td(:,2)=cellstr([string(idx); string(idx); string(idx)]);
        end

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

        t.Operations = [reduced_operations];
        %t.Operations = [reduced_operations reduced_operations reduced_operations];

        save(strcat(BASE_OUTPUT_FOLDER, filesep, STAGE_SAVE_FILENAME), '-struct', 't');

        %% Load timeseries information
        ts = load(hctsafile,'TimeSeries');
        ts = ts.TimeSeries;
        tst = struct2table(ts);

        segment_idx=[];
        segment_cluster=[];
        for i = 1:no_of_channels
            x = length(idx);
            segment_idx = [segment_idx, [(i-1)*x+1:i*x]];
            segment_cluster = [segment_cluster, idx'];
        end

        %%
        segment_cluster_summary=array2table([segment_idx',(segment_idx')*epoch_seconds,segment_cluster']);

        set(gcf,'Visible','on')
        load(EDF_FILE);
        unique_clusters = unique(segment_cluster);
        for c = 1:size(unique_clusters,2)
            cluster_label = unique_clusters(c);
            seg_idxs = find(find(table2array(segment_cluster_summary(:, 3)) == cluster_label) <= size(datamat, 1));
            sub_clusterdatamat = datamat(seg_idxs, :);

            [s, s_optk, ch, ch_optk, db, db_optk] = calculate_opt_k(sub_clusterdatamat);
            opt_k_report_str = strcat(opt_k_report_str, ...
                sprintf("%s\n------------\n(Number of epochs=%d) \nSilhouette: %.04f (%d) \nCalinskiHarabasz: %.04f (%d) \nDaviesBouldin: %.04f (%d)\n\n", ...
            sprintf("Cluster %d:", cluster_label), ...
            size(sub_clusterdatamat, 1), ...
            s, ...
            s_optk, ...
            ch, ...
            ch_optk, ...
            db, ...
            db_optk));
    %         [idx, c, sse] = perform_clustering(sub_clusterdatamat, s_opt, epoch_seconds);
    %         plot_cluster_fft(header, data, hctsafile, epoch_seconds, timeseries_sampling_rate, seg_idxs, ...
    %             idx, sprintf('Cluster %s', str2double(cluster_label))); 
        end
        %%

        %% Hierarchical clustering
        % eucD=pdist(datamat, 'cosine');
        % clustTreeEuc=linkage(eucD, 'average');
        % 
        % disp(cophenet(clustTreeEuc, eucD));
        % 
        % [h,nodes] = dendrogram(clustTreeEuc,100);
        % h_gca = gca;
        % h_gca.TickDir = 'out';
        % h_gca.TickLength = [.002 0];
        % h_gca.XTickLabel = [];

        % [h,nodes] = dendrogram(clustTreeEuc,20);
        % hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',.185);
        % idx=hidx;
        % n_clust = length(unique(hidx));

        set(gcf,'Visible','on')

        clearvars -except SBJ_ID opt_k_report_str BASE_OUTPUT_FOLDER config;
    end;

    fid = fopen(strcat(BASE_OUTPUT_FOLDER, filesep, "bipolar_1_EEG_Main_opt_k.txt"), 'wt');
    fprintf(fid, opt_k_report_str);
    fclose(fid);

end
%% Functions

function [idx, c, sse] = perform_clustering(datamat, n_clust, epoch_seconds)
    [idx, c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);
                    
    segment_idx=[];
    segment_cluster=[];
    for i = 1:1
        x = length(idx);
        segment_idx = [segment_idx, [(i-1)*x+1:i*x]];
        segment_cluster = [segment_cluster, idx'];
    end

    %%
    segment_cluster_summary=array2table([segment_idx',(segment_idx')*epoch_seconds,segment_cluster']);
    %writetable(segment_cluster_summary, [output_folder filesep 'cluster_segment_information.csv']);
    
    unique_clusters = unique(segment_cluster);
    for c = 1:size(unique_clusters,2)
        cluster_label = unique_clusters(c);
        seg_idxs = find(table2array(segment_cluster_summary(:, 3)) == cluster_label);
        sub_clusterdatamat = datamat(seg_idxs, :);
        
        [s, s_optk, ch, ch_optk, db, db_optk] = calculate_opt_k(sub_clusterdatamat);
    end    
end


