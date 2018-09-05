%% Configuration
% dataset - col 1 = data set number, col 2 = start epoch, col 3 = end epoch
datasets = [ ...
    1, 334, 1374; ...
    5, 380, 1442; ...
    7, 391, 1207; ...
   13, 375, 1495; ...
   14, 174, 1353; ...
   36, 436, 1435; ...
   68, 302, 1316; ... % Epoch 30, 31, 32 were removed. (Not used because mismatch with annotation)
   78, 243, 1307; ...
  129, 160, 1311; ...
  159, 345, 1388; ...
  180, 353, 1440; ...
  224, 228, 1266; ...
  246,  91, 1050; ...
  294, 199, 1243; ...
  302, 124, 1329; ...
  439, 150, 1164; ...
  458, 277, 1374; ...
  557, 328, 1421; ... % Not used because all data are mostly NaN.
  579, 194, 1258; ...
  596, 266, 1396; ...
  604, 502, 1457; ...
  658, 291, 1354; ...
  684, 243, 1271; ...
  687, 332, 1365; ...
  748, 488, 1191; ...
  749,  69, 1009; ...
  752, 147, 1096; ...
  774, 183, 1216; ...
  807, 329, 1232; ...
  813, 286, 1179; ...
  821, 314, 1316; ...
  869, 343, 1405; ...
  870, 282, 1286; ...
];

DATASET=870;
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
TARGET_FOLDER=strcat('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800', num2str(DATASET), '_HCTSA_Analysis');
TARGET_FILE=strcat('HCTSA_1800', num2str(DATASET), '_N.mat');
ANSWER_FILE=strcat('/Volumes/Spaceship/ccshs_datasets/ccshs_1800', num2str(DATASET), '_annot.mat');
NUM_OF_CHANNELS=7;
TARGET_HCTSA_FILE='HCTSA_EEG_EOG_N.mat';
MODE=1; % Mode = 0 (Split file), Mode = 1 (Ran specific plot)
DRAW_PLOT=1;

% Derived
datasets=datasets';
validData = datasets(1, :); % Valid datasets
if ~ismember(DATASET,validData)
    error(strcat('Dataset is invalid...input: ', validData, '\n'));
end
endW = datasets(2, :);
endS = datasets(3, :);
endID = find(DATASET==validData); % Determine which ending/beginning to use

START = endW(endID)+1;
END = endS(endID)-1;

%% Main execution body
homedir=pwd;
%% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)

%% Read HCTSA top 200 features
hctsafile=strcat(TARGET_FOLDER, filesep, TARGET_FILE);
all_op = load(strcat(TARGET_FOLDER, filesep, TARGET_FILE),'Operations');

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
    
    if (match_input_index > 0)
        nn=nn+1;
        feat_id(nn) = match_input_index;
        feat(nn).id = match_input_index; % all_op.Operations(i).ID % Actual operation ID
        feat(nn).name = op_name;
    end
end
clear i n nn op_name name

%% Duplicate the original file
cd(TARGET_FOLDER)
copyfile(TARGET_FILE, 'HCTSA_N.mat');

annotation = load(ANSWER_FILE);
sleepstage = annotation.sleepstage(START:END);

%% Mode 0 - Split the HCTSA based on channels and extract only the top 200 features
if (MODE == 0)
    norm_struct = load('HCTSA_N.mat');
    norm_ts = struct2table(norm_struct.TimeSeries);
    splt_tbl=split(string(table2array(norm_ts(:,1))), "_");
    unique_channels=unique(splt_tbl(:,1), 'stable');
    no_of_epochs_per_channels=size(norm_ts,1)/size(unique_channels, 1);

    for i=1:size(unique_channels, 1)
        chn_struct = norm_struct;
        chn_struct.TS_DataMat = norm_struct.TS_DataMat(no_of_epochs_per_channels*(i-1)+START:no_of_epochs_per_channels*(i-1)+END, feat_id);
        chn_struct.TS_Quality = norm_struct.TS_Quality(no_of_epochs_per_channels*(i-1)+START:no_of_epochs_per_channels*(i-1)+END, feat_id);
        chn_struct.TimeSeries = norm_struct.TimeSeries(no_of_epochs_per_channels*(i-1)+START:no_of_epochs_per_channels*(i-1)+END);

        chn_ts_table = struct2table(chn_struct.TimeSeries);
        chn_ts_table(:,2) = array2table(cellstr(string(sleepstage)));
        chn_ts_table(:,5) = array2table([START:END]');
        chn_struct.TimeSeries = table2struct(chn_ts_table);

        chn_struct.Operations = norm_struct.Operations(feat_id);
        chn_struct.op_clust.ord = 1:size(feat_id, 2);
        chn_struct.ts_clust.ord = 1:size(chn_struct.TS_DataMat, 1);
        
        % Save to default HCTSA.mat so that annotation can be labelled.
        save(strcat('HCTSA.mat'), '-struct', 'chn_struct');
        TS_LabelGroups([]);

        movefile('HCTSA.mat', char(strcat('HCTSA_', unique_channels(i), '_N.mat')));
    end
    
    %% Now we construct the combination of channels
    % EEG + EOG
    combo_channels={[1, 2], [1, 3], [1, 2,3]};
    combo_channel_names={"EEG_EOG", "EEG_EMG", "EEG_EOG_EMG"};
    for i=1:length(combo_channels)
        combo_channel=combo_channels{i};
        combo_name=combo_channel_names{i};

        ts_indexes = [];
        ss = [];
        for j=1:length(combo_channel)
            ii = combo_channel(j);
            ts_indexes = [ts_indexes, no_of_epochs_per_channels*(ii-1)+START:no_of_epochs_per_channels*(ii-1)+END];
            ss = [ss; sleepstage]; 
        end
        
        chn_struct = norm_struct;
        chn_struct.TS_DataMat = norm_struct.TS_DataMat(ts_indexes, feat_id);
        chn_struct.TS_Quality = norm_struct.TS_Quality(ts_indexes, feat_id);
        chn_struct.TimeSeries = norm_struct.TimeSeries(ts_indexes);
        
        chn_ts_table = struct2table(chn_struct.TimeSeries);
        chn_ts_table(:,2) = array2table(cellstr(string(ss)));
        chn_struct.TimeSeries = table2struct(chn_ts_table);

        chn_struct.Operations = norm_struct.Operations(feat_id);
        chn_struct.op_clust.ord = 1:(size(feat_id, 2) * length(combo_channel));
        chn_struct.ts_clust.ord = 1:size(chn_struct.TS_DataMat, 1);
        
        % Save to default HCTSA.mat so that annotation can be labelled.
        save(strcat('HCTSA.mat'), '-struct', 'chn_struct');
        TS_LabelGroups([]);

        movefile('HCTSA.mat', char(strcat('HCTSA_', combo_name, '_N.mat')));        
    end
    
end

%% Mode = 1 (Ran plot)
if (MODE == 1)
    disp(strcat('Processing-', TARGET_HCTSA_FILE, '...'));
    
    copyfile(TARGET_HCTSA_FILE, 'HCTSA_N.mat');

    t=load('HCTSA_N.mat');
    ts=struct2table(t.TimeSeries);
    ts_labels = table2array(ts(:,1))';
    
    FIGURE_PREFIX=extractBefore(TARGET_HCTSA_FILE, ".mat");

    %%
    if (DRAW_PLOT == 1)
        %%
        % PLOT 1: Plot that preserve the ordering of timeseries and features
        TS_plot_DataMatrix('norm');
        dcm=datacursormode(gcf);
        dcm.UpdateFcn = {@custom_cursor_text, sleepstage, feat, ts_labels};
        savefig(strcat(FIGURE_PREFIX, '_Unclustered.fig'));

        wake_idx = find(sleepstage==0);
        n1_idx = find(sleepstage==1);
        n2_idx = find(sleepstage==2);
        n3_idx = find(sleepstage==3);
        rem_idx = find(sleepstage==5);

        reorder_sleepstages = [zeros(1, length(wake_idx)),  ones(1, length(n1_idx)), ...
            ones(1, length(n2_idx)) * 2, ...
            ones(1, length(n3_idx)) * 3, ...
            ones(1, length(rem_idx)) * 5];
        reorder_based_on_annotation = [wake_idx; n1_idx; n2_idx; n3_idx; rem_idx];
        TS_plot_DataMatrix('norm', 'customOrder', {reorder_based_on_annotation', []});
        dcm=datacursormode(gcf);
        dcm.UpdateFcn = {@custom_cursor_text, reorder_sleepstages, feat, ts_labels(reorder_based_on_annotation)};

        sleep_stages_count = [length(wake_idx), length(n1_idx), length(n2_idx), length(n3_idx), length(rem_idx)];
        draw_lines(feat_id, sleep_stages_count);
        savefig(strcat(FIGURE_PREFIX, '_Unclustered_With_SleepStage.fig'));

        %%
        % Cluster the data so that we can display group information
        distanceMetricRow = 'euclidean';
        linkageMethodRow = 'average';
        distanceMetricCol = 'corr_fast';
        linkageMethodCol = 'average';

        TS_cluster(distanceMetricRow, linkageMethodRow, distanceMetricCol, linkageMethodCol);

        %%
        % PLOT 2: Cluster information (After TS_cluster)
        TS_plot_DataMatrix('cl', 'colorGroups', false);
        dcm=datacursormode(gcf);
        dcm.UpdateFcn = {@custom_cursor_text, sleepstage, feat, ts_labels};
        savefig(strcat(FIGURE_PREFIX, '_Clustered_NoGroup.fig'));

        TS_plot_DataMatrix('cl', 'colorGroups', false, 'customOrder', {reorder_based_on_annotation', []});
        dcm=datacursormode(gcf);
        dcm.UpdateFcn = {@custom_cursor_text, reorder_sleepstages, feat, ts_labels(reorder_based_on_annotation)};
        draw_lines(feat_id, sleep_stages_count);
        savefig(strcat(FIGURE_PREFIX, '_Clustered_NoGroup_With_SleepStage.fig'));

        %%
        % PLOT 3: Group information (After TS_cluster)
        TS_plot_DataMatrix('cl', 'colorGroups', true);
        dcm=datacursormode(gcf);
        dcm.UpdateFcn = {@custom_cursor_text, sleepstage, feat, ts_labels};
        savefig(strcat(FIGURE_PREFIX, '_Clustered_Group.fig'));

        TS_plot_DataMatrix('cl', 'colorGroups', true, 'customOrder', {reorder_based_on_annotation', []});
        dcm=datacursormode(gcf);
        dcm.UpdateFcn = {@custom_cursor_text, reorder_sleepstages, feat, ts_labels(reorder_based_on_annotation)};
        draw_lines(feat_id, sleep_stages_count);
        savefig(strcat(FIGURE_PREFIX, '_Clustered_Group_With_SleepStage.fig'));
    end
    
    %%
    % TS_Classify
    TS_classify('norm', 'svm_linear', 1);
    
    %% 
    % Top_Features
    copyfile('HCTSA_N.mat', 'HCTSA.mat');
    TS_TopFeatures();
    
    
    
end
    
%% 
cd(homedir);

%% Function definitions

function draw_lines(feat_id, sleep_stages_count)
    fs = size(feat_id, 2);
    line([0, fs], ones(1, 2) .* sum(sleep_stages_count(1:1)), 'LineWidth', 1, 'Color', 'black');
    text(fs + 1, sum(sleep_stages_count(1:1)),'W');
    
    line([0, fs], ones(1, 2) .* sum(sleep_stages_count(1:2)), 'LineWidth', 1, 'Color', 'black');
    text(fs + 1, sum(sleep_stages_count(1:2)),'N1');

    line([0, fs], ones(1, 2) .* sum(sleep_stages_count(1:3)), 'LineWidth', 1, 'Color', 'black');
    text(fs + 1, sum(sleep_stages_count(1:3)),'N2');

    line([0, fs], ones(1, 2) .* sum(sleep_stages_count(1:4)), 'LineWidth', 1, 'Color', 'black');
    text(fs + 1, sum(sleep_stages_count(1:4)),'N3');
    
    line([0, fs], ones(1, 2) .* sum(sleep_stages_count(1:5)), 'LineWidth', 1, 'Color', 'black');
    text(fs + 1, sum(sleep_stages_count(1:5)),'R');
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
