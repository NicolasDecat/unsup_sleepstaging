%% Configuration

SBJ_ID='KJ_N1';
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
TARGET_FOLDER=strcat('/Volumes/Spaceship/Voss_Lucid/', SBJ_ID, '/ALL_EEG');
NUM_SECONDS=30;

MODE='MAIN_CLUSTER';
%MODE='SUB_CLUSTER';
%MODE='PLOT_MAIN_CLUSTER';
%MODE='PLOT_SUB_CLUSTER';

PLOT_SUB_CLUSTER_MAX_SAMPLE=10;
PLOT_RANDOM=1;
PLOT_ONLY_CLUSTER=[1];

TARGET_FILE=strcat('HCTSA_N_', SBJ_ID, '_1_EEG_Main_5_Clusters.mat');

% SUB_TARGET_FILE='HCTSA_N_1_EEG_6_REM_substages.mat';
% SUB_TARGET_FILE='HCTSA_N_1_EEG_6_REM_second_level_substages.mat';
% SUB_TARGET_FILE='HCTSA_N_ME_N2_1_EEG_5_REM_substages.mat';
% SUB_TARGET_FILE='HCTSA_N_ME_N2_1_EEG_Main.mat';
SUB_TARGET_FILE = strcat('HCTSA_N_', SBJ_ID, '_Cluster_4_1_EEG_2_REM_substages.mat');
%SUB_TARGET_FILE = strcat('HCTSA_N_', SBJ_ID, '_1_EEG_7_REM_second_level_substages.mat');

BASE_PLOT_FOLDER=strcat(TARGET_FOLDER, filesep, extractBefore(SUB_TARGET_FILE, '.'));
if strcmp(MODE, 'PLOT_MAIN_CLUSTER')
    BASE_PLOT_FOLDER=strcat(TARGET_FOLDER, filesep, extractBefore(TARGET_FILE, '.'));
end

EDF_FILE=strcat('/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/', SBJ_ID, '_data.mat');

%potential_lucid_dreaming_epoch_start = [17550, 20940, 21210, 22200, 22440, 22530];
potential_lucid_dreaming_epoch_start = [];

%% Main execution body
distanceMetricRow = 'euclidean'; %
linkageMethodRow = 'average'; %
distanceMetricCol = 'corr_fast';
linkageMethodCol = 'average'; %

homedir=pwd;
%% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)

%% Features
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

%%
cd(TARGET_FOLDER)

warning('off');
if (strcmp(MODE, 'MAIN_CLUSTER'))
    % Copy the current file to standard HCTSA.mat naming for LabelGroups
    copyfile(TARGET_FILE, 'HCTSA.mat');

    TS_LabelGroups([]);

    copyfile('HCTSA.mat', 'HCTSA_N.mat');

    %% Sort based on algorithm clusters
    t=load('HCTSA_N.mat');
    ts=struct2table(t.TimeSeries);
    datamat=t.TS_DataMat;
    ts_labels = table2array(ts(:,1))';
    ts_algo_groups = str2num(char(table2array(ts(:,2))));

    g1_idx = find(ts_algo_groups==1);
    g2_idx = find(ts_algo_groups==2);
    g3_idx = find(ts_algo_groups==3);
    g4_idx = find(ts_algo_groups==4);
    g5_idx = find(ts_algo_groups==5);
    reorder_sleepstages = [ones(1, length(g1_idx)),  ones(1, length(g2_idx)) * 2, ...
        ones(1, length(g3_idx)) * 3, ...
        ones(1, length(g4_idx)) * 4, ...
        ones(1, length(g5_idx)) * 5];

    %%
    % Cross reference with another sample data (sleep.org) to calculate the
    % correlation between the features
    %
    %
    % REF_DATA = '/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800749_HCTSA_Analysis/HCTSA_C3-A2_N.mat';
    % ref=load(REF_DATA);
    % rs=struct2table(ref.TimeSeries);
    % datamat_ref=ref.TS_DataMat;
    % ts_ref_groups = str2num(char(table2array(rs(:,2))));
    % ts_ref_groups=ts_ref_groups+1;
    % ts_ref_groups(ts_ref_groups==6) = 5;
    % 
    % corm=zeros(5,5);
    % for i = 1:5
    %     datamat_single = datamat(find(ts_algo_groups==i), :);
    %     
    %     corr_value = [];
    %     for j = 1:5
    %         datamat_ref_single = datamat_ref(find(ts_ref_groups==j), :);
    %         c=corrcoef(mean(datamat_single), mean(datamat_ref_single(:, 1:195)));
    %         corm(i, j) = c(1,2);
    %     end
    % end
    % 
    % disp(corm);
    % [a, b] = max(abs(corm'));
    % disp(a);
    % disp(b);

    %%
    reorder_based_on_annotation = [g1_idx; g2_idx; g3_idx; g4_idx; g5_idx];
    sleep_stages_count = [length(g1_idx), length(g2_idx), length(g3_idx), length(g4_idx), length(g5_idx)];
    reorder_ts_labels = ts_labels(reorder_based_on_annotation);

%     TS_plot_DataMatrix('norm', 'customOrder', {reorder_based_on_annotation', []});
%     draw_lines(feat_id, sleep_stages_count);
%     dcm=datacursormode(gcf);
%     dcm.UpdateFcn = {@custom_cursor_text, reorder_sleepstages, feat, reorder_ts_labels};
%     plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, reorder_ts_labels);

    %% Perform HCTSA clustering    
    TS_cluster(distanceMetricRow, linkageMethodRow, distanceMetricCol, linkageMethodCol);

    %%
    % HCTSA clusters
    tsl = TS_plot_DataMatrix('cl', 'colorGroups', false);
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;
    %plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, tsl');

    %%
    tsl = TS_plot_DataMatrix('cl', 'colorGroups', true);
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;
    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, tsl');

    %%
    hctsa_cluster_labels=tsl;
    %hctsa_cluster_labels=split(string(tsl), '_');
    %hctsa_cluster_labels = double(hctsa_cluster_split(:, 3));
    
    hctsa_reorder_based_on_annotation = [];
    for k=1:length(hctsa_cluster_labels)
        hctsa_reorder_based_on_annotation = [hctsa_reorder_based_on_annotation, ...
            find(hctsa_cluster_labels==string(reorder_ts_labels(k)))];
    end
    
    tl=TS_plot_DataMatrix('cl', 'colorGroups', false, 'customOrder', {hctsa_reorder_based_on_annotation', []});
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;

    draw_lines(feat_id, sleep_stages_count);
    dcm=datacursormode(gcf);
    dcm.UpdateFcn = {@custom_cursor_text, reorder_sleepstages, feat, reorder_ts_labels};
    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, reorder_ts_labels);

    %%
    TS_plot_DataMatrix('cl', 'colorGroups', true, 'customOrder', {hctsa_reorder_based_on_annotation', []});
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;

    draw_lines(feat_id, sleep_stages_count);
    dcm=datacursormode(gcf);
    dcm.UpdateFcn = {@custom_cursor_text, reorder_sleepstages, feat, reorder_ts_labels};
    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, reorder_ts_labels);


elseif strcmp(MODE, 'SUB_CLUSTER')
    
    copyfile(SUB_TARGET_FILE, 'HCTSA.mat');
    TS_LabelGroups([]);

    copyfile('HCTSA.mat', 'HCTSA_N.mat');

    %% Sort based on algorithm clusters
    t=load('HCTSA_N.mat');
    ts=struct2table(t.TimeSeries);
    datamat=t.TS_DataMat;
    ts_labels = table2array(ts(:,1))';
    ts_algo_groups = str2num(char(table2array(ts(:,2))));
    
    uniq_values = unique(ts_algo_groups);
    
    reorder_substages = [];
    reorder_based_on_annotation = [];
    sleep_stages_count = [];
    
    for i=1:length(uniq_values)
        g_idx = find(ts_algo_groups==uniq_values(i));
        reorder_substages = [reorder_substages, ones(1, length(g_idx)) * i];
        reorder_based_on_annotation = [reorder_based_on_annotation; g_idx];
        sleep_stages_count = [sleep_stages_count; length(g_idx)];
    end

    reorder_ts_labels = ts_labels(reorder_based_on_annotation);

    %% Perform HCTSA clustering
    TS_cluster(distanceMetricRow, linkageMethodRow, distanceMetricCol, linkageMethodCol);

    %%
    TS_plot_DataMatrix('norm', 'customOrder', {reorder_based_on_annotation', []});
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;

    draw_lines(feat_id, sleep_stages_count);
    dcm=datacursormode(gcf);
    dcm.UpdateFcn = {@custom_cursor_text, reorder_substages, feat, reorder_ts_labels};
    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, reorder_ts_labels);

    %%
    % HCTSA clusters
    tsl = TS_plot_DataMatrix('cl', 'colorGroups', false);
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;

    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, tsl');

    tsl = TS_plot_DataMatrix('cl', 'colorGroups', true);
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;

    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, tsl');

    %%
    hctsa_cluster_split=split(string(tsl), '_');
    hctsa_cluster_labels = double(hctsa_cluster_split(:, 3));
    
    hctsa_reorder_based_on_annotation = [];
    ts_indexes = [];
    for k=1:length(hctsa_cluster_labels)
        ts_idx = find(endsWith(table2array(ts(:,1)), strcat('_', num2str(hctsa_cluster_labels(k)))));
        ts_indexes = [ts_indexes, ts_idx];
    end
    
    for k=1:length(ts_indexes)
        hctsa_reorder_based_on_annotation = [hctsa_reorder_based_on_annotation, ...
            find(ts_indexes==reorder_based_on_annotation(k))];
    end
    
    TS_plot_DataMatrix('cl', 'colorGroups', false, 'customOrder', {hctsa_reorder_based_on_annotation', []});
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;

    draw_lines(feat_id, sleep_stages_count);
    dcm=datacursormode(gcf);
    dcm.UpdateFcn = {@custom_cursor_text, reorder_substages, feat, reorder_ts_labels};
    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, reorder_ts_labels);

    %%
    TS_plot_DataMatrix('cl', 'colorGroups', true, 'customOrder', {hctsa_reorder_based_on_annotation', []});
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;
    
    draw_lines(feat_id, sleep_stages_count);
    dcm=datacursormode(gcf);
    dcm.UpdateFcn = {@custom_cursor_text, reorder_substages, feat, reorder_ts_labels};
    plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, reorder_ts_labels);    

elseif strcmp(MODE, 'PLOT_SUB_CLUSTER') || strcmp(MODE, 'PLOT_MAIN_CLUSTER') 

    load(EDF_FILE);

    %% Extract channel names, sampling frequency from header struct
    channel = cellstr(header.channelname);  % Channel names + converted into cell array
    fs = header.samplerate(1);      % Assume sampling rate are the same for all channels
    recordtime = header.records;    % aaTotal recorded time (seconds)
    nchannels = header.channels;    % Number of channels

    % Amplitude calibration
    phys_range = header.physmax - header.physmin;
    dig_range = header.digimax - header.digimin;

    gain = phys_range ./ dig_range;

    for i = 1:size(data, 1)
       data(i, :) = (data(i, :) - header.digimin(i)) * gain(i) + header.physmin(i);
    end

    channel_info=[(1:length(header.channelname))',string(header.channelname)]

    cd(homedir);
    kmeans_clustering_configuration_channels;
    cd(TARGET_FOLDER);

    %% Pre-processing

    % Notch-filter of 50Hz

    % Design a filter with a Q-factor of Q=35 to remove a 50 Hz tone from 
    % system running at 300 Hz.
    Wo = 50/(fs/2);  BW = Wo/35;
    [b,a] = iirnotch(Wo,BW);
    data(chans.EEG,:)=filter(b,a,data(chans.EEG,:));

    % Bandpass between [0.5 70]
    [b,a]=butter(2,[0.5,70]/(fs/2),'bandpass');
    sub_ts=data(chans.EEG,:)';
    sub_ts=filtfilt(b,a,sub_ts);
    data(chans.EEG,:) = sub_ts';

    [b,a]=butter(2,[0.5,15]/(fs/2),'bandpass');
    sub_ts=data(chans.EOG,:)';
    sub_ts=filtfilt(b,a,sub_ts);
    data(chans.EOG,:) = sub_ts'; 
    
    ds=[];
    for i=1:size(data,1)
       ds=[ds; downsample(data(i, :), 2)];
    end
    data=ds;
    % data = medfilt1(data,3);
    fs = fs/2;
    es = NUM_SECONDS;
    
    %% Plot the relevant channel
    HCTSA_to_use = SUB_TARGET_FILE;
    if strcmp(MODE, 'PLOT_MAIN_CLUSTER') 
        HCTSA_to_use = TARGET_FILE;
    end
    
    t=load(HCTSA_to_use);
    dm = t.TS_DataMat;
    quality = t.TS_Quality;
    ts = t.TimeSeries;
    original_ts = struct2table(t.TimeSeries);
    single_channel_size = size(dm,1)/6;

    clusters = table2array(unique(original_ts(:,2)))';
    %%
    % Output CSV
    for i = 1:length(clusters)
        table_csv = original_ts;
        name_split = split(string(table2array(table_csv(:, 1))), '_');
        table_csv(:, 1) = cellstr(name_split(:, 3));
    end
    
    %%
    table_csv_trim = table_csv(:, [1, 2]);
    extra_column = array2table((str2double(table2array(table_csv_trim(:,1)))-1)*30, 'VariableNames', {'EpochStarts'});
    table_csv_trim = [table_csv_trim extra_column];
    writetable(table_csv_trim, strcat(TARGET_FOLDER, filesep, SUB_TARGET_FILE, '.csv'));
    
    %% Plot substage epochs
    set(gcf,'Visible','off')
    for i = 1:length(clusters)
        cluster_foldername = strcat(BASE_PLOT_FOLDER, filesep, 'Cluster', clusters(i));
        mkdir(char(cluster_foldername));
        
        if (length(PLOT_ONLY_CLUSTER) ~= 0 && isempty(find(PLOT_ONLY_CLUSTER == str2double(clusters(i)))))
           continue; 
        end
        
        cluster_ts = original_ts(original_ts.Keywords == string(clusters(i)), :);
            
        set(gcf,'Visible','off');

        if (PLOT_SUB_CLUSTER_MAX_SAMPLE ~= 0 && PLOT_RANDOM)
            max_rows = size(cluster_ts, 1);
            if (max_rows >= PLOT_SUB_CLUSTER_MAX_SAMPLE)
                max_rows = PLOT_SUB_CLUSTER_MAX_SAMPLE;
            end
            rand_ts = randperm(size(cluster_ts, 1), max_rows);
            cluster_ts = cluster_ts(rand_ts, :);
        end
        
        number_count = 0;
        for j = 1:size(cluster_ts, 1)

            if PLOT_SUB_CLUSTER_MAX_SAMPLE ~= 0 && number_count >= PLOT_SUB_CLUSTER_MAX_SAMPLE
                break;
            end
            
            f=figure('rend','painters','pos',[10 10 1200 800]);
            set(f,'color','w');
            
            clust_label = cluster_ts(j, 1);
            clust_label_split = split(string(table2array(clust_label)), '_');
            
            tseg = str2double(clust_label_split(3));
            startIndex = (str2double(clust_label_split(3)) - 1) * NUM_SECONDS;
            endIndex = startIndex + NUM_SECONDS;
            
            draw(f, channel_info, data, startIndex, fs, es, chans, sprintf('EEG, EOG and EMG between %d and %d seconds', startIndex, endIndex));

            set(gcf,'Visible','off');

            imagename = strcat('LD_Clust_Sub_', num2str(tseg,'%04d'),'.png');
            saveas(gcf,string(strcat(cluster_foldername, filesep, imagename)));
            close
            
            number_count = number_count + 1;
        end
    end
end    
%%
cd(homedir);

function draw(f, channel_info, data, num_seconds, fs, es, chans, t)
    no_of_epochs=(num_seconds/es)+1;
    EEG_CHN = chans.EEG;
    EOG_CHN = chans.EOG;
    EMG_CHN = chans.EMG;

    ax=subplot(3, 1, 1);
    draw_subplot(ax, data, "EEG", EEG_CHN, cellstr(channel_info(EEG_CHN, 2))', no_of_epochs, [-200 200], fs, es, chans);
    title(t);
    
    ax=subplot(3, 1, 2);
    draw_subplot(ax, data, "EOG", EOG_CHN, cellstr(channel_info(EOG_CHN, 2))', no_of_epochs, [-50, 50], fs, es, chans);

    ax=subplot(3, 1, 3);
    draw_subplot(ax, data, "EMG", EMG_CHN, cellstr(channel_info(EMG_CHN, 2))', no_of_epochs, [-10, 10], fs, es, chans);
end

function draw_subplot(ax, data, plot_title, chnls, channel_labels, no_of_epochs, y_limits, fs, es, chans)
    pos=get(ax, 'position');
    pos(4) = 0.26;
    set(ax, 'position', pos);
    
    for i=1:length(chnls)
        chn_idx=chnls(i);
        ts_data = data(chn_idx,:);

        frame=es*fs;
        
        if chn_idx==15 || chn_idx==16
            plt=plot(ax, frame*(no_of_epochs-1)+1:frame*no_of_epochs, ts_data(:, frame*(no_of_epochs-1)+1:frame*no_of_epochs), 'LineWidth', 0.5);
        else
            plt=plot(ax, frame*(no_of_epochs-1)+1:frame*no_of_epochs, ts_data(:, frame*(no_of_epochs-1)+1:frame*no_of_epochs), 'LineWidth', 0.5);
        end
        hold on;

        axis tight;
        ax.YGrid = 'on';
        ax.XGrid = 'on';

        pos = get(ax, 'Position');
        pos(1) = 0.055;
        pos(3) = 0.9;
        set(ax, 'Position', pos)

        set(ax,'box','off');
        ylabel(ax, plot_title);
        legend(channel_labels);
        legend(ax, 'show');
        set(ax, 'XTick', 0:fs:size(data, 2));
        XTickLabel = get(ax,'XTick');
        set(ax, 'XTickLabel', string(XTickLabel/fs));
        
%         YTickLabel = get(ax,'YTick');
%         set(ax, 'YTickLabel', string(YTickLabel));

        if (length(y_limits) > 0)
            %ax.YAxisLocation = 'origin';
            ylim(ax, y_limits);
        end
    end
end


function draw_potential_lucid_lines(feat_id, yaxis)
    fs = size(feat_id, 2);
    line([0, fs], ones(1, 2) .* yaxis, 'LineWidth', 1, 'Color', 'red');
    text(fs + 1, yaxis,'LD', 'Color', 'red');
end

function plot_lucid_dreaming_candidates(potential_lucid_dreaming_epoch_start, feat_id, reorder_ts_labels)
    if length(potential_lucid_dreaming_epoch_start) == 0
        return;
    end
    
    for i = 1:length(potential_lucid_dreaming_epoch_start)
       time_seg = (potential_lucid_dreaming_epoch_start(i)+NUM_SECONDS)/NUM_SECONDS;
       epoch_no_idx = find(endsWith(string(reorder_ts_labels), strcat('_timeseg_',num2str(time_seg)))==1);
       if (~isempty(epoch_no_idx))
        draw_potential_lucid_lines(feat_id, epoch_no_idx);
       end
    end
end

function draw_lines(feat_id, sleep_stages_count)
    fs = size(feat_id, 2);
    
    for i=1:length(sleep_stages_count)
        line([0, fs], ones(1, 2) .* sum(sleep_stages_count(1:i)), 'LineWidth', 1, 'Color', 'black');
        text(fs + 1, sum(sleep_stages_count(1:i)),strcat('C', num2str(i)));
    end
end


function output_txt = custom_cursor_text(obj, event_obj, ss, fea, tl)
    % Display the position of the data cursor
    % obj          Currently not used (empty)
    % event_obj    Handle to event object
    % output_txt   Data cursor text (character vector or cell array of character vectors).

    pos = get(event_obj,'Position');
    stage = ss(pos(2));
    feature = fea(pos(1)).name;
    ts_label = char(string(tl(pos(2))));
    
    output_txt = {['X: ',num2str(pos(1),4)],...
        ['Y: ',num2str(pos(2),4)], ...
        ['Time Series: ', ts_label], ...
        ['Stage: ', num2str(stage)], ...
        ['Feature: ', feature], ...
        };

    % If there is a Z-coordinate in the position, display it as well
    if length(pos) > 2
        output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
    end
end