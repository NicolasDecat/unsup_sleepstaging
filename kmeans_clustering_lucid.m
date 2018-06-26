clear feat_id feat_id_test;

% %% CONFIGURATION

% Global
n_cluster_samples=5;
n_clust = 5;
n_channels = 1;

all_3_eegs=0;

%% Configuration for testing (only applicable if perform_test=1)
perform_test=0;
test_hctsafile='/Volumes/Spaceship/Voss_Lucid/KJ_N1/30seconds/HCTSA_N.mat';
% TEST_C4_COL=1:789;
% TEST_EOG_COL=(length(TEST_C4_COL)*3+1):((length(TEST_C4_COL)*3+1)+length(TEST_C4_COL)-1);
% TEST_EMG_COL=(length(TEST_C4_COL)*5+1):((length(TEST_C4_COL)*5+1)+length(TEST_C4_COL)-1);

% test_hctsafile='TS_SK_N1_F1_30sec_6chan/HCTSA_N.mat';
% TEST_C4_COL=1:20;
% TEST_EOG_COL=(length(TEST_C4_COL)*3+1):((length(TEST_C4_COL)*3+1)+length(TEST_C4_COL)-1);
% TEST_EMG_COL=(length(TEST_C4_COL)*5+1):((length(TEST_C4_COL)*5+1)+length(TEST_C4_COL)-1);
% TEST_K = 10;

kmeans_clustering_configuration;
%C4_COL=F4_COL;

%%
all_op = load(hctsafile,'Operations');
all_op_test = load(test_hctsafile, 'Operations');

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
    match_input_index = 0;
    match_test_index = 0;
    
    for i = 1:length(all_op.Operations)
        name = all_op.Operations(i).Name;
        if strcmp(op_name,name)
            match_input_index = i;
            break;
        end
    end
    
    for i = 1:length(all_op_test.Operations)
        name = all_op_test.Operations(i).Name;
        if strcmp(op_name,name)
            match_test_index = i;
            break;
        end
    end
    
    % If the feature was found in both data set, then include them
    if perform_test
        if (match_input_index > 0 & match_test_index > 0)
            nn=nn+1;
            feat_id(nn) = match_input_index;
            feat(nn).id = match_input_index; % all_op.Operations(i).ID % Actual operation ID
            feat(nn).name = op_name;
            
            feat_id_test(nn) = match_test_index;
            feat_test(nn).id = match_test_index; % all_op.Operations(i).ID % Actual operation ID
            feat_test(nn).name = op_name;
        end
    else
        if (match_input_index > 0)
            nn=nn+1;
            feat_id(nn) = match_input_index;
            feat(nn).id = match_input_index; % all_op.Operations(i).ID % Actual operation ID
            feat(nn).name = op_name;
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

%% Find the max cluster that has the highest silhouette values
% max_cluster = 0;
% max_s = 0;
% for i = 1:20
%     [idx,c, sse] = kmeans(datamat,i,'Distance','sqeuclidean',...
%                         'Display','off','Replicates',50,'MaxIter',500);
%     s=silhouette(datamat, idx);
%     if (mean(s) > max_s)
%         max_s = mean(s);
%         max_cluster = i;
%     end
% end
% 
% disp(['Max cluster is ' num2str(max_cluster) ' with ' num2str(max_s)]);

%% Perform k-means clustering
% for i=5:30
%     n_clust=i;
    [idx, c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);
                    
    potential_lucid_dreaming_stages=idx(592:598);
    disp(potential_lucid_dreaming_stages);
    
    
    %% Save HTCSA_N trimmed for HCTSA analysis (assigned labels)
    t=load('/Volumes/Spaceship/Voss_Lucid/KJ_N1/HCTSA/HCTSA_N.mat');
    dm = t.TS_DataMat;
    quality = t.TS_Quality;
    ts = t.TimeSeries;

    single_channel_size = size(dm,1)/6;

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
    t.Operations = [reduced_operations reduced_operations reduced_operations];

    save('/Volumes/Spaceship/Voss_Lucid/KJ_N1/ALL_EEG/HCTSA_N_1_EEG_TEST.mat', '-struct', 't');

    
    %% Lucid dreaming substage
    % Supposely the most count of potential lucid dreaming is the REM
    % stage. Use that and split the REM stage into two clusters.
    rem_stage = mode(potential_lucid_dreaming_stages);
    idx_rem_stage = find(idx==rem_stage);
    
    n_clust = 2;
    datamat=datamat(idx_rem_stage,:);
    [idx, c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);    
    
    disp(idx(find(ismember(idx_rem_stage, (592:598)))));
   
    
    %% Save HTCSA_N trimmed (substages)
    t=load('/Volumes/Spaceship/Voss_Lucid/KJ_N1/HCTSA/HCTSA_N.mat');
    dm = t.TS_DataMat;
    quality = t.TS_Quality;
    ts = t.TimeSeries;
    original_ts = struct2table(t.TimeSeries);

    single_channel_size = size(dm,1)/6;

    startIndex=1; % No trim

    % 1 EEG channel
    t.TS_DataMat = [dm(idx_rem_stage, feat_id)];
    t.TS_Quality = [quality(idx_rem_stage, feat_id)];
    t.TimeSeries = [ts(idx_rem_stage)];

    ts = t.TimeSeries;
    td=struct2table(ts);
    td(:,2)=cellstr([string(idx)]);
    ts=table2struct(td);
    
    t.TimeSeries = ts;
    
    selected_timeseries = original_ts(C4_COL, :);
    EOG_timeseries = original_ts(EOG_COL,:);
    EMG_timeseries = original_ts(EMG_COL, :);
    
    %% Plot substage epochs
    for i = 1:n_clust
        selected_eeg_cluster = selected_timeseries(idx==i,:);
        selected_eog_cluster = EOG_timeseries(idx==i,:);
        selected_emg_cluster = EMG_timeseries(idx==i,:);
        random_cluster_ids = table2array(selected_eeg_cluster(:, 1));

        for j = 1:length(random_cluster_ids)
            random_cluster_id = random_cluster_ids(j);
            index = find(table2array(selected_eeg_cluster(:, 1)) == string(random_cluster_id));

            single_signal=cell2mat(table2array(selected_eeg_cluster(index, 4)));
            signals{1}.signal=single_signal;
            signals{1}.signal_type='EEG';
            signals{1}.min = -1700;
            signals{1}.max = 1700;
            single_signal_length=length(single_signal);

            single_signal=cell2mat(table2array(selected_eog_cluster(index, 4)));
            signals{2}.signal=single_signal;
            signals{2}.signal_type='EOG';
            signals{2}.min = -3950;
            signals{2}.max = 3950;

            single_signal=cell2mat(table2array(selected_emg_cluster(index, 4)));
            signals{3}.signal=single_signal;
            signals{3}.signal_type='EMG';
            signals{3}.min = 600;
            signals{3}.max = -600;

            tslabel = string(table2cell(selected_eeg_cluster(index, 2)));
            labels = split(tslabel, ',');
            timeseg = labels(2);
            timeseg = strrep(timeseg, 'timeseg_', '');
            tseg = str2double(timeseg);
            timeseg_seconds_start = (tseg-1)*epoch_seconds;    
            timeseg_seconds_end = tseg*epoch_seconds;

            set(gcf,'Visible','off')
            draw_figure(epoch_seconds, signals, timeseries_sampling_rate, epoch_seconds, (size(ts,1)/no_of_channels)*single_signal_length, ...
                sprintf('Cluster: %d TS: %s Time: between %d secs and %d secs' , i, timeseg, timeseg_seconds_start, timeseg_seconds_end));

            set(gcf,'Visible','off');

            imagename = strcat('LD_Clust_',num2str(i,'%01d'),'_',num2str(j),'_',num2str(tseg,'%04d'),'.png');
            saveas(gcf,[output_folder filesep imagename]) % saveas, imwrite or imsave? print(imagename,'-dpng')?
            close
        end
    end


    %%

    % 3 EEG channels
%     t.TS_DataMat = [dm(idx_rem_stage, feat_id) dm(idx_rem_stage+single_channel_size, feat_id) dm(idx_rem_stage+(single_channel_size*2), feat_id)];
%     t.TS_Quality = [quality(idx_rem_stage, feat_id) quality(idx_rem_stage+single_channel_size, feat_id) quality(idx_rem_stage+(single_channel_size*2), feat_id)];
%     t.TimeSeries = [ts(idx_rem_stage);ts(idx_rem_stage+single_channel_size);ts(idx_rem_stage+(single_channel_size*2))];
% 
%     ts = t.TimeSeries;
%     td=struct2table(ts);
%     td(:,2)=cellstr([string(idx); string(idx); string(idx)]);
%     ts=table2struct(td);
%     
%     t.TimeSeries = ts;

    % 2 channels
    % t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(startIndex+single_channel_size:single_channel_size*2, feat_id)];
    % t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(startIndex+single_channel_size:single_channel_size*2, feat_id)];
    % t.TimeSeries = [ts(startIndex:single_channel_size);ts(startIndex+single_channel_size:single_channel_size*2)];

    % 3 channels
    % t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(single_channel_size+startIndex:single_channel_size*2, feat_id) dm(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
    % t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(single_channel_size+startIndex:single_channel_size*2, feat_id) quality(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
    % t.TimeSeries = [ts(startIndex:single_channel_size);ts(single_channel_size+startIndex:single_channel_size*2);ts(single_channel_size*2+startIndex:single_channel_size*3)];

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

%     t.Operations = [reduced_operations reduced_operations reduced_operations];
    t.Operations = [reduced_operations];
    save('/Volumes/Spaceship/Voss_Lucid/KJ_N1/ALL_EEG/HCTSA_N_1_EEG_2_substage_TEST.mat', '-struct', 't');

    
%%     
%     if (length(unique(idx(592:598))) == 1 & length(unique(idx(716:735))) ==  1 & ...
%             idx(591)~=idx(592) & idx(598) ~= idx(599) & ...
%             idx(715)~=idx(716) & idx(735) ~= idx(736))
%         disp(['Cluster ' num2str(i) ' has same consecutive clusters for lucid']);
%     end
% end

% [s,~]=silhouette(datamat, idx);
% disp(mean(s));
% n_clust=max_cluster;

%% Hierarchical clustering
% eucD=pdist(datamat, 'cosine');
% clustTreeEuc=linkage(eucD, 'average');
% 
% disp(cophenet(clustTreeEuc, eucD));
% 
% [h,nodes] = dendrogram(clustTreeEuc,0);
% h_gca = gca;
% h_gca.TickDir = 'out';
% h_gca.TickLength = [.002 0];
% h_gca.XTickLabel = [];
% 
% [h,nodes] = dendrogram(clustTreeEuc,20);
% hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',.185);
% idx=hidx;
% n_clust = length(unique(hidx));

%% Plot

% markersize=12;         
% feature1 = 1;
% feature2 = 2;
% 
set(gcf,'Visible','off');

% Generate random colors
% colors = distinguishable_colors(n_clust);
% figure;
% 
% for i=1:n_clust
%     plot(datamat(idx==i,feature1),datamat(idx==i,feature2),'.','color', colors(i,:),'MarkerSize',markersize);
%     legend_str{i} = ['Cluster ' num2str(i) '(SSE: ' num2str(sse(i)) ')'];
%     hold on;
% end
% 
% dcm_obj1 = datacursormode(gcf);
% set(dcm_obj1,'UpdateFcn',{@cluster_plot_data_cursor,datamat(:, feature1:feature2),idx});
% 
% title(sprintf("KMeans Clustering - Lucid - Chan: %s", "C4"));
% plot(c(:,feature1),c(:,feature2),'kx','MarkerSize',markersize,'LineWidth',2);
% %plot(c(:,1),c(:,2),'ko','MarkerSize',12,'LineWidth',2)
% 
% legend1 = legend(string(legend_str), 'Location','NW');
% set(legend1,'Location','northeastoutside');   
% 
% set(gcf,'Visible','off')
% 
% savefig([output_folder filesep 'sp_feature1_feature2.fig']);

%% Load timeseries information
ts = load(hctsafile,'TimeSeries');
ts = ts.TimeSeries;
tst = struct2table(ts);
selected_timeseries = tst(C4_COL, :);
eeg_mat = cell2mat(table2array(selected_timeseries(:,4)));

EOG_timeseries = tst(EOG_COL,:);
eog_mat = cell2mat(table2array(selected_timeseries(:,4)));

EMG_timeseries = tst(EMG_COL, :);
emg_mat = cell2mat(table2array(selected_timeseries(:,4)));

segment_idx=[];
segment_cluster=[];
for i = 1:1
    x = length(idx);
    segment_idx = [segment_idx, [(i-1)*x+1:i*x]];
    segment_cluster = [segment_cluster, idx'];
end

segment_cluster_summary=array2table([segment_idx',(segment_idx')*epoch_seconds,segment_cluster']);
writetable(segment_cluster_summary, [output_folder filesep 'cluster_segment_information.csv']);

%% Perform validation test (using another dataset).
if perform_test
    
    datamat_test = load(test_hctsafile,'TS_DataMat');
    orig_datamat_test = datamat_test.TS_DataMat;

    if all_3_eegs
        TEST_O2_COL=TEST_C4_COL+length(TEST_C4_COL);
        TEST_F4_COL=TEST_O2_COL+length(TEST_C4_COL);
        datamat_test = [orig_datamat(TEST_C4_COL,:) orig_datamat(TEST_O2_COL,:) orig_datamat(TEST_F4_COL,:)];
    else
        datamat_test = orig_datamat_test(TEST_C4_COL,:);
    end

    if n_channels == 1
        datamat_test = datamat_test(:,feat_id_test);
    elseif n_channels == 2
        datamat_test = [datamat_test(:,feat_id_test) orig_datamat_test(TEST_EOG_COL,feat_id_test)];
    elseif n_channels == 3
        datamat_test = [datamat_test(:,feat_id_test) orig_datamat_test(TEST_EOG_COL,feat_id_test) orig_datamat_test(TEST_EMG_COL,feat_id_test)];        
    end
    
    for n=1:size(datamat_test,1)
        % Calculate Euclidean distance from datapoint to each centre
        distance = sqrt(sum((bsxfun(@minus,c,datamat_test(n,:))).^2,2));
        % Find the cluster of minimum distance
        [~, idx_test(n)] = min(distance);

        % Using nearest neighbours
%         test_idx = knnsearch(datamat, datamat_test(n,:), 'k', TEST_K)';
%         idx_test(n) = mode(idx(test_idx));
    end
    
    test_cluster_summary=array2table([(1:size(datamat_test,1))', ((1:size(datamat_test,1))*epoch_seconds)', idx_test']);
    writetable(test_cluster_summary, [output_folder filesep 'test_segment_information.csv']);
end

if n_cluster_samples > 0
    for i = 1:n_clust
        selected_eeg_cluster = selected_timeseries(idx==i,:);
        selected_eog_cluster = EOG_timeseries(idx==i,:);
        selected_emg_cluster = EMG_timeseries(idx==i,:);
%         random_cluster_ts_idx = randperm(size(selected_eeg_cluster, 1));
% 
%         if (length(random_cluster_ts_idx) <= n_cluster_samples)
%             random_cluster_ids = random_cluster_ts_idx(1:length(random_cluster_ts_idx));
%         else
%             random_cluster_ids = random_cluster_ts_idx(1:n_cluster_samples);
%         end
        random_cluster_ids = table2array(selected_eeg_cluster(:, 1));

        for j = 1:length(random_cluster_ids)
            random_cluster_id = random_cluster_ids(j);
            index = find(table2array(selected_eeg_cluster(:, 1)) == string(random_cluster_id));

            single_signal=cell2mat(table2array(selected_eeg_cluster(index, 4)));
            signals{1}.signal=single_signal;
            signals{1}.signal_type='EEG';
            signals{1}.min = -1700;
            signals{1}.max = 1700;
            single_signal_length=length(single_signal);

            single_signal=cell2mat(table2array(selected_eog_cluster(index, 4)));
            signals{2}.signal=single_signal;
            signals{2}.signal_type='EOG';
            signals{2}.min = -3950;
            signals{2}.max = 3950;

            single_signal=cell2mat(table2array(selected_emg_cluster(index, 4)));
            signals{3}.signal=single_signal;
            signals{3}.signal_type='EMG';
            signals{3}.min = 600;
            signals{3}.max = -600;

            tslabel = string(table2cell(selected_eeg_cluster(index, 2)));
            labels = split(tslabel, ',');
            timeseg = labels(2);
            timeseg = strrep(timeseg, 'timeseg_', '');
            tseg = str2double(timeseg);
            timeseg_seconds_start = (tseg-1)*epoch_seconds;    
            timeseg_seconds_end = tseg*epoch_seconds;

            set(gcf,'Visible','off')
            draw_figure(epoch_seconds, signals, timeseries_sampling_rate, epoch_seconds, (size(ts,1)/no_of_channels)*single_signal_length, ...
                sprintf('Cluster: %d TS: %s Time: between %d secs and %d secs' , i, timeseg, timeseg_seconds_start, timeseg_seconds_end));

            set(gcf,'Visible','off');

            imagename = strcat('LD_Clust_',num2str(i,'%01d'),'_',num2str(j),'_',num2str(tseg,'%04d'),'.png');
            saveas(gcf,[output_folder filesep imagename]) % saveas, imwrite or imsave? print(imagename,'-dpng')?
            close
        end
    end
end

set(gcf,'Visible','on')

% keywords_col = [keywords_col, array2table(idx)];
% timescale=(mod(table2array(keywords_col(:,2)), 120)*5);
% timescale(timescale==0) = 600;
% keywords_col = [keywords_col, array2table(timescale)];

% draw_figure(30, selectedSignal, 200);

function draw_figure(tmax, selectedSignal, samplingRate, epoch_seconds, whole_ts_length, title_array)
    % Get number of signals
    num_signals = length(selectedSignal);
    timeID=1;
    x0=0;
    y0=0;
    width=1800;
    height=900;
    
    figure;
    for s = 1:num_signals
        % Get signal
        signal =  selectedSignal{s}.signal;
%         signal = downsample(signal, 2);
%         samplingRate = samplingRate/2;
        
        t = [0:whole_ts_length-1]/samplingRate; % = record_duration

        % Parameters for normalisation - use global max and min if amplitude
        % matters. If not, set an arbitary value
%         sigMin = -0.3; %min(signal);
%         sigMax = 0.3; %max(signal);
        sigMin = selectedSignal{s}.min;
        sigMax = selectedSignal{s}.max;
        signalRange = sigMax - sigMin;
    %     
        % Identify indexes of 30 seconds of signal according to tstart, tend
        % Otherwise, indexes = find(t<=tmax);
        tStart = find(t==(timeID-1)*epoch_seconds);
        tEnd = find(t==timeID*epoch_seconds)-1;
        indexes = tStart:1:tEnd;
        signal = signal(indexes);
        time = t(1:length(indexes)); % time = t(indexes); % Hide real time, always display 0 -30 seconds 

    %% Switch-case for num_signals
         if signalRange~= 0
            signal = signal/(sigMax-sigMin);
         end
         
    switch (num_signals)
        case 1
            % Centred around 0
            signal = signal - mean(signal);
        case {2,3}
            % Add signal below the previous one
            signal = signal - mean(signal) + (num_signals - s + 1);
            %     signal = signal + (num_signals - s + 1); % Without zero-centred
            %     signal = signal - 0.5*mean(signal) + (num_signals - s + 1);
            % Plot line dividing signals
            plot(time,s-0.5*ones(1,length(time)),'color',[0.5,0.5,0.5])
    end

    % Color code signal type - can be customised + will depends on screen setting
    switch (selectedSignal{s}.signal_type)
        case 'EEG'
            ccode = [0.1,0.5,0.8];
        case 'EOG'
            ccode = [0.1,0.5,0.3];
        case 'EMG'
            ccode = [0.8,0.5,0.2];
        case 'ECG'
            ccode = [0.8,0.1,0.2];
        otherwise
            ccode = [0.2,0.2,0.2];
    end
    set(gcf,'units','pixels','position',[x0,y0,width,height]);
    plot(time, signal,'Color',ccode);
    hold on;
end
    
    %% Plot configuration
    grid on
    ax = gca;
    fig = gcf;
    
    switch (num_signals)
        case 1
            % Set axes limits
            v = axis();
            v(1:2) = [0,tmax];
            v(3:4) = [sigMin, sigMax];
            
            axis(v);
            % Set x-axis 
            xlabel('Time(sec)')
            ax.XTick = [0:epoch_seconds];
            ax.FontSize = 10;
            % Set y-axis
            ylabel('Amplitude(\muV)')
            ax.YTick = linspace(sigMin,sigMax,epoch_seconds);

        case {2,3}
            % Set axis limits
            v = axis();
            v(1:2) = [0,tmax];
            v(3:4) = [0.5 num_signals+0.5];
            axis(v);
            % Set x axis
            xlabel('Time(sec)');
            ax = gca;
            ax.XTick = [0:epoch_seconds];
            ax.FontSize = 10;
            % Set y axis labels
            ylabel('Amplitude (mV)')
            %% Without scale
            signalLabels = cell(1,num_signals); %Revert the order such that first channel stays on top
            % for s = 1:num_signals
            %     signalLabels{num_signals-s+1} = selectedHeader(s).signal_labels;
            % end
            %ax.YTick = 1:num_signals;
            %ax.YTickLabels = signalLabels;
            %% With scale
            ax.YTick = [0.55,1,1.45,1.55,2,2.45,2.55,3,3.45];
            ax.YTickLabels = {'-300',selectedSignal{3}.signal_type,'+300','-300',selectedSignal{2}.signal_type,'+300','-300',selectedSignal{1}.signal_type,'+300'};
            ax.FontSize = 15;
            % 
            % Set figure size
            fig.Units = 'pixels';
            fig.Position = [x0,y0,width, height];
            fig.Color = [0.95 0.95 0.95];
            title(title_array);

    end

    % Reduce white space ** Can be adjusted
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + 1.1*ti(1);
    bottom = outerpos(2) + 1.1*ti(2);
    ax_width = outerpos(3) - 1.1*ti(1) - 3*ti(3);
    ax_height = outerpos(4) - 1.1*ti(2) - 4*ti(4);
    ax.Position = [left bottom ax_width ax_height];

end

function txt = cluster_plot_data_cursor(~,event_obj,myarray, indexes)
    pos = get(event_obj,'Position');
    
    idx = find((abs(myarray(:,1)-pos(1))<0.0001) & (abs(myarray(:,2)-pos(2))<0.0001));
    timeseconds = (idx-1)*30;
    txt = {['X: ',num2str(pos(1))],...
           ['Y: ',num2str(pos(2))],...
           ['Index: ',num2str(idx)],...
           ['Cluster: ', num2str(indexes(idx))],...
           ['Time: ' , num2str(timeseconds)]};
end



