clear feat_id feat_id_test;

%n_clust = 3; % Uncomment this to FIX cluster.
for n_clust = 2:6

% %% CONFIGURATION

% Global
n_cluster_samples=1;
substage_samples=1;
n_channels = 1;

all_3_eegs=0;

SUB_ID='ME_N3';
SUB_SECONDARY_ID='_bipolar';

%% Configuration for testing (only applicable if perform_test=1)
kmeans_clustering_configuration;

STAGE_SAVE_FILENAME=strcat('HCTSA_N_', SUB_ID, SUB_SECONDARY_ID, '_1_EEG_Main_', num2str(n_clust), '_Clusters.mat');
EDF_FILE=strcat('/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/', SUB_ID, '_data.mat');
%BASE_OUTPUT_FOLDER='/Volumes/Spaceship/Voss_Lucid/KJ_N1/ALL_EEG';
BASE_OUTPUT_FOLDER=strcat('/Volumes/Spaceship/Voss_Lucid/', SUB_ID, SUB_SECONDARY_ID, '/1_CHANS');
%%
all_op = load(hctsafile,'Operations');
OPS_FILE='reduced_ops.txt';
%%
currentdir = pwd;

fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);

%% Wanted operation names
feat_name = features{1,2};

%% Check operation name, get feat_id
nn=0;
for n = 1:length(feat_name)
    op_name = char(feat_name(n));
    for i = 1:length(all_op.Operations)
        name = all_op.Operations(i).Name;
        if strcmp(op_name,name)
            nn=nn+1;
            feat_id(nn) = i;
            feat(nn).id = i; % all_op.Operations(i).ID % Actual operation ID
            feat(nn).name = cellstr(name);
        end
    end
end
clear i n nn op_name name

%% Use feat_id to select data from full op
datamat = load(hctsafile,'TS_DataMat');
orig_datamat = datamat.TS_DataMat;

if all_3_eegs
    O2_COL=C4_COL+length(C4_COL);
    F4_COL=O2_COL+length(C4_COL);
    datamat = [orig_datamat(C4_COL, feat_id) orig_datamat(O2_COL, feat_id) orig_datamat(F4_COL, feat_id)];
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

calculate_opt_k(datamat, "Main cluster: ");

%% Perform k-means clustering
% for i=5:30
%     n_clust=i;
    [idx, c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);
                    
%     potential_lucid_dreaming_stages=idx(592:598);
%     disp(potential_lucid_dreaming_stages);
    
    %% Save HTCSA_N trimmed for HCTSA analysis (assigned labels)
    t=load(hctsafile);
    dm = t.TS_DataMat;
    quality = t.TS_Quality;
    ts = t.TimeSeries;

    single_channel_size = size(dm,1)/no_of_channels;

    startIndex=1; % No trim
    endIndex=single_channel_size;

    mkdir(output_folder);

    %1 EEG channel
    t.TS_DataMat = [dm(startIndex:endIndex, feat_id)];
    t.TS_Quality = [quality(startIndex:endIndex, feat_id)];
    t.TimeSeries = [ts(startIndex:endIndex)];

    %3 EEG channels
%     t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(single_channel_size+startIndex:single_channel_size*2, feat_id) dm(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
%     t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(single_channel_size+startIndex:single_channel_size*2, feat_id) quality(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
%     t.TimeSeries = [ts(startIndex:single_channel_size);ts(single_channel_size+startIndex:single_channel_size*2);ts(single_channel_size*2+startIndex:single_channel_size*3)];

    %2 channels
%     t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(startIndex+single_channel_size:single_channel_size*2, feat_id)];
%     t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(startIndex+single_channel_size:single_channel_size*2, feat_id)];
%     t.TimeSeries = [ts(startIndex:single_channel_size);ts(startIndex+single_channel_size:single_channel_size*2)];

    %3 channels
%     t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(single_channel_size+startIndex:single_channel_size*2, feat_id) dm(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
%     t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(single_channel_size+startIndex:single_channel_size*2, feat_id) quality(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
%     t.TimeSeries = [ts(startIndex:single_channel_size);ts(single_channel_size+startIndex:single_channel_size*2);ts(single_channel_size*2+startIndex:single_channel_size*3)];

    ts = t.TimeSeries;
    td=struct2table(ts);
    %td(:,2)=cellstr([string(idx); string(idx); string(idx)]);
    %td(:,2)=cellstr([string(idx); string(idx); string(idx); string(idx); string(idx); string(idx)]);
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

    t.Operations = [reduced_operations];
    %t.Operations = [reduced_operations reduced_operations reduced_operations];

    save(strcat(BASE_OUTPUT_FOLDER, filesep, STAGE_SAVE_FILENAME), '-struct', 't');

    %% Load timeseries information
    ts = load(hctsafile,'TimeSeries');
    ts = ts.TimeSeries;
    tst = struct2table(ts);

    segment_idx=[];
    segment_cluster=[];
    for i = 1:1
        x = length(idx);
        segment_idx = [segment_idx, [(i-1)*x+1:i*x]];
        segment_cluster = [segment_cluster, idx'];
    end

    %%
    segment_cluster_summary=array2table([segment_idx',(segment_idx')*epoch_seconds,segment_cluster']);
    writetable(segment_cluster_summary, [output_folder filesep 'cluster_segment_information.csv']);
    
    set(gcf,'Visible','on')
    load(EDF_FILE);
    unique_clusters = unique(segment_cluster);
    for c = 1:size(unique_clusters,2)
        cluster_label = unique_clusters(c);
        seg_idxs = find(table2array(segment_cluster_summary(:, 3)) == cluster_label);
        sub_clusterdatamat = datamat(seg_idxs, :);
        
        [s, s_opt, ~, ~, ~, ~] = calculate_opt_k(sub_clusterdatamat, sprintf("Cluster %d:", cluster_label));
        
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

clear all;
end;
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


