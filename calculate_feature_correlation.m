% Select features from reduced_ops.txt, make another .mat file to work on
% cross-validation.
%configuration_settings;

% HCTSA_FILE='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/HCTSA_N.mat';
%HCTSA_FILE='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800001_HCTSA/1800001_EOG_5chan/HCTSA_N.mat';

all_dataset = [ ...
  % Male = 1 (First 5) Male = 0 (Rest)
  "1800439", 197;...
  "1800749", 196;...
  "1800752", 194;...
  "1800807", 197;...
  "1800870", 198;...
  "1800458", 197;...
  "1800596", 197;...
  "1800604", 198;...
  "1800748", 197;...  
  "1800821", 195;...
 ];

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

datasets = datasets';   

number_of_channels=7;
xtra_correlation_channels = [...
  "C3-A2", "LOC-A2"; ...
  "C3-A2", "EMG1-EMG2"; ...
  "LOC-A2", "EMG1-EMG2"; ...
  "LOC-A2", "ROC-A1"; ...
  "LOC-A2+ROC-A1", "LOC-A2-ROC-A1"; ...
  "LOC-A2+ROC-A1", "LOC-A2*ROC-A1"; ...
  "LOC-A2-ROC-A1", "LOC-A2*ROC-A1"; ...
];

table_columns = {'Dataset', 'Channel', 'Mean_r', 'Mean_abs_r', 'Mean_r2'};
table = array2table(zeros(0,length(table_columns)));
table.Properties.VariableNames = table_columns;

eog_grand_corr_matrix = zeros(4);

for dataset_count=1:size(all_dataset, 1)
    dataset_name = all_dataset(dataset_count, 1);
    dataset_feature_size = str2num(char(all_dataset(dataset_count, 2)));
    
    DATASET=str2num(char(strrep(dataset_name, '1800', '')));
    DATASET_NAME=char(dataset_name);
    submatrix_size=dataset_feature_size;

    HCTSA_FILE=strcat('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_', DATASET_NAME, '_HCTSA/HCTSA_N.mat');
    OUTPUT_FILE=strcat(DATASET_NAME, '_feature_corr.mat');
    OPS_FILE='reduced_ops.txt';

    validData = datasets(1, :); % Valid datasets
    if ~ismember(DATASET,validData)
        % whichData = 1; % Default case 
        % *** can change to display error and keep asking for valid input using
        % Add while loop?
        error(strcat('Dataset is invalid...input: ', validData, '\n'));
    end
    %% Remove initial W stage from randomisation and sampling
    % Marking the end of W stage and Sleep stages (N1-N3) (Manually)
    endW = datasets(2, :);
    endS = datasets(3, :);
    endID = find(DATASET==validData); % Determine which ending/beginning to use

    START = endW(endID)+1;
    END = endS(endID)-1;

    %addpath('/usr/local/gurobi/7.5.1/matlab');

    % homedir = pwd;
    % %% Start HCTSA tools
    % cd(HCTSA_DIR)
    % startup
    % cd(homedir)

    %% Read text file
    fileID = fopen(OPS_FILE);
    features = textscan(fileID,'%s %s %s');
    fclose(fileID);

    %% Wanted operation names
    feat_name = features{1,2};

    %% All operation names
    hctsafile = HCTSA_FILE;
    all_op = load(hctsafile,'Operations');

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
    h = load(hctsafile);

    datamat = h.TS_DataMat;

    [timeseries,features]=size(datamat);
    hctsa_ops = datamat(:,feat_id);

    single_channel_size=size(hctsa_ops,1)/number_of_channels;

    t=struct2table(h.TimeSeries);
    s=split(table2array(t(:,1)), '_');
    c=unique(s(:,1), 'stable')

    % CHANGEME
    % eeg_ops=hctsa_ops(1:single_channel_size,:);
    % eog_ops=hctsa_ops(single_channel_size+1:single_channel_size*2,:);
    % emg_ops=hctsa_ops(single_channel_size*2+1:single_channel_size*3,:);

    channels={};

    % leog_ops=hctsa_ops(1:single_channel_size,:);
    % reog_ops=hctsa_ops(single_channel_size+1:single_channel_size*2,:);
    % l_plus_r_ops=hctsa_ops(single_channel_size*2+1:single_channel_size*3,:);
    % l_minus_r_ops=hctsa_ops(single_channel_size*3+1:single_channel_size*4,:);
    % l_times_r_ops=hctsa_ops(single_channel_size*4+1:single_channel_size*5,:);

    for i = 1:number_of_channels
        a_channel = hctsa_ops(single_channel_size*(i-1)+1:single_channel_size*i,:);
        channels{i} = a_channel;
    end

    all_channels_ops=[];

    for i = 1:number_of_channels
        channels{i} = channels{i}(START:END, :);
        all_channels_ops = [all_channels_ops channels{i}];
    end

    %%
    % Combine all features
    %all_channels_ops = [eeg_ops eog_ops emg_ops];

    % SELECT_TOP_200_FEATURES=[8,14,19,22,25,26,27,29,39,41,45,61,73,77,94,95,97,98,99,121,122,123,138,144,150,153,154,155,159,164,172,210,218,223,228,240,243,254,255,267,268,270,273,276,277,281,286,293,296,297,311,314,318,319,320,322,330,341,342,343,345,357,358,367,368,370,372,384,388,408,415,416,443,444,461,488,572];
    % 
    % all_channels_ops(:, SELECT_TOP_200_FEATURES);

    %% Iterate through channels
    for chn = 1:number_of_channels
        channel_name=string(c(chn));
        PLOT_LABEL_PREFIX=channel_name;
        PLOT_TITLE=channel_name;

        [corrl, pv]=corrcoef(channels{chn});
        corr_r2=corrl.^2;

        %%
        A = corr_r2;

        Y=A;
        Y(logical(eye(size(Y)))) = 0; 

        m = [mean(mean(corrl)), mean(mean(abs(corrl))), mean(mean(corr_r2))];
        row = [DATASET_NAME, channel_name, m(1), m(2), m(3)];
        table = [table; array2table(row, 'VariableNames', table_columns)];
        
%         fprintf('\nMean r: %.4f\n', m(1));
%         fprintf('\nMean abs(r): %.4f\n', m(2));
%         fprintf('\nMean r^2: %.4f\n', m(3));
% 

        %% GUROBI
        % CM=Y;
        % num_to_keep=submatrix_size;
        % 
        % % corr=corrcoef(emg_ops);
        % n=size(corr, 1);
        % varnames=cellstr([strcat('EEG-', string(feat_id)) strcat('EOG-', string(feat_id)) strcat('EMG-', string(feat_id))]);
        % imagesc(corr); % plot the matrix
        % set(gca, 'XTick', 1:n); % center x-axis ticks on bins
        % set(gca, 'YTick', 1:n); % center y-axis ticks on bins
        % set(gca, 'XTickLabel', varnames); % set x-axis labels
        % set(gca, 'YTickLabel', varnames); % set y-axis labels
        % title(sprintf('Correlation Matrix between %dx%d feature matrix for EEGxEOGxEMG (Dataset 1)',n, n), 'FontSize', 14); % set title
        % colormap('jet'); % set the colorscheme
        % colorbar; % enable colorbar

        %% 
        clear corr;

        % eeg_values=find_top_X_features(eeg_ops, 90);
        % eog_values=find_top_X_features(eeg_ops, 100)+198;
        % emg_values=find_top_X_features(eeg_ops, 10)+(198*2);
        % 
        % for i=1:length(eeg_values)
        %    fprintf('%d,',eeg_values(i)); 
        % end
        % 
        % for i=1:length(eog_values)
        %    fprintf('%d,',eog_values(i)); 
        % end
        % 
        % for i=1:length(emg_values)
        %    fprintf('%d,',emg_values(i)); 
        % end


        %%
%         figure;
%         n=size(corrl, 1);
%         varnames=cellstr([strcat(PLOT_LABEL_PREFIX, '-', string(feat_id))]);
%         imagesc(corrl); % plot the matrix
%         set(gca, 'XTick', 1:n); % center x-axis ticks on bins
%         set(gca, 'YTick', 1:n); % center y-axis ticks on bins
%         set(gca, 'XTickLabel', varnames); % set x-axis labels
%         set(gca, 'YTickLabel', varnames); % set y-axis labels
%         set(gcf, 'Units', 'Inches', 'Position', [0, 0, 10.125, 8.125], 'PaperUnits', 'Inches', 'PaperSize', [10.125, 8.125]);
%         % set(gcf, 'Units', 'Inches', 'Position', [0, 0, 12.125, 10.125], 'PaperUnits', 'Inches', 'PaperSize', [12.125, 10.125]);
%         xtickangle(90);
%         title(sprintf('Correlation Matrix between %dx%d feature matrix for %s (Dataset %d)',n, n, PLOT_TITLE, DATASET), 'FontSize', 19); % set title
%         colormap('jet'); % set the colorscheme
%         colorbar; % enable colorbar

        %     N = size(CM,1); 
        % 
        %     clear model;  
        %     names = strseq('x',[1:N]);  
        %     model.varnames = names;  
        %     model.Q = sparse(CM); % Gurobi needs a sparse matrix as input  
        %     model.A = sparse(ones(1,N));  
        %     
        %     obj = zeros(1,N);
        %     model.obj = obj;
        %     model.rhs = num_to_keep;  
        %     model.sense = '=';  
        %     model.vtype = 'B';
        % 
        %     gurobi_write(model, 'qp.mps');
        % 
        %     results = gurobi(model);
        %     
        % 
        %     values = results.x;
        %     save(OUTPUT_FILE, 'corr_eeg', 'corr_eog', 'corr_emg', 'submatrix_size', 'model', 'results', 'values');
    end
    
    %% Calculate cross-channels correlations
    for xtra_channel = 1:size(xtra_correlation_channels, 1)
        channel_1 = xtra_correlation_channels(xtra_channel, 1);
        channel_2 = xtra_correlation_channels(xtra_channel, 2);

        first_chn_idx = find(string(c)==channel_1);
        second_chn_idx = find(string(c)==channel_2);
        [corrm, p]=corrcoef(channels{first_chn_idx}, channels{second_chn_idx});
        corrl = corrm(2, 1);
        corr_r2=corrl.^2;
        m = [corrl abs(corrl) corr_r2];
        row = [DATASET_NAME, strcat(channel_1, 'x', channel_2), m(1), m(2), m(3)];
        table = [table; array2table(row, 'VariableNames', table_columns)];
    end
    
    %% Calculate cross-channels correlations (single features)
    for xtra_channel = 1:size(xtra_correlation_channels, 1)
        channel_1 = xtra_correlation_channels(xtra_channel, 1);
        channel_2 = xtra_correlation_channels(xtra_channel, 2);

        first_chn_idx = find(string(c)==channel_1);
        second_chn_idx = find(string(c)==channel_2);
        
        feature_corr = [];
        chn_1 = channels{first_chn_idx};
        chn_2 = channels{second_chn_idx};
        
        for f = 1:size(channels{first_chn_idx}, 2)
            [corrm, p] = corrcoef(chn_1(:,f), chn_2(:, f));
            corrl = corrm(2, 1);
            corr_r2=corrl.^2;
            feature_corr = [feature_corr; [corrl, abs(corrl), corr_r2]];
        end
        
        m = [mean(feature_corr(:,1)) mean(feature_corr(:,2)) mean(feature_corr(:,3))];
        row = [DATASET_NAME, strcat(channel_1, '#', channel_2), m(1), m(2), m(3)];
        table = [table; array2table(row, 'VariableNames', table_columns)];
    end    
    
    %% Calculate the relationship between dfferent channels per feature
    eeg=channels{1};
    leog=channels{2};
    emg=channels{3};
    reog=channels{4};
    plus=channels{5};
    minus=channels{6};
    times=channels{7};

    cm = zeros(4);
    cm2 = zeros(4);
    
    for f = 1:size(channels{first_chn_idx}, 2)
        t = [leog(:,f) reog(:,f) minus(:, f) times(:, f)];
        cm = cm + corrcoef(t);
        cm2 = cm2 + (corrcoef(t).^2);
    end
    
    eog_grand_corr_matrix = eog_grand_corr_matrix + (cm./size(channels{first_chn_idx}, 2));
%     cm2./size(channels{first_chn_idx}, 2)
end

%table

%% Calculate the mean of all correlations
[UA, ~, idx] = unique(table(:,[2]), 'stable');
accum_mean_corr = [UA, ...
    array2table(accumarray(idx,double(table2array(table(:,3))),[],@mean), 'VariableNames', { 'Mean_r' }) ...
    array2table(accumarray(idx,double(table2array(table(:,4))),[],@mean), 'VariableNames', { 'Mean_abs_r' }) ...
    array2table(accumarray(idx,double(table2array(table(:,5))),[],@mean), 'VariableNames', { 'Mean_r2' }) ...
    ]

%% Construct error bar
standard_errors = [];
mean_of_errors = [];
for i = 1:size(UA, 1)
    t = table(table.Channel == table2array(UA(i, 1)), :);
    r2 = str2double(t.Mean_r2);
    standard_errors = [standard_errors, std(r2)./sqrt(size(all_dataset, 1))];
    mean_of_errors = [mean_of_errors; mean(r2)];
end

%% Plot 1 - Correlation between features
feature_index=find(~contains(table2array(accum_mean_corr(:,1)), 'x') & ~contains(table2array(accum_mean_corr(:,1)), '#'));
accum_table=accum_mean_corr(feature_index,:);
unique_type = unique(accum_table(:, 1), 'stable');
figure('pos',[250 250 1200 800]);

bar_data = [];
for i = 1:size(unique_type, 1)
   filter = accum_table(accum_table.Channel == table2array(unique_type(i, :)), :);

   bar_data = [bar_data; table2array(filter(:, 4))'];
end

bar_data = bar_data;
standard_errors = standard_errors;

%standard_error_reshape = reshape(standard_errors, 3, 6)';
barplot = barwitherr(standard_errors(feature_index)', bar_data);

% barplot = bar(bar_data);
% l = cell(1,3);
% l{1} = 'EEG';
% l{2} = 'EEG + EOG';
% l{3} = 'EEG + EOG + EMG';
% legend(barplot, l);

title(strcat('Mean Correlation of determination (R^2) between features in each channel (n=', num2str(size(all_dataset, 1)), ')'), 'FontSize', 16, 'Interpreter', 'tex');
grid on;
xlabel('Channel', 'FontSize', 15);
text(1:length(bar_data),bar_data,num2str(bar_data),'vert','bottom','horiz','center', 'FontSize', 14); 
ylabel('Mean Correlation of determination (R^2)', 'FontSize', 15, 'Interpreter', 'tex');
set(gca, 'xticklabel', table2array(unique_type));
%fix_xticklabels(gca, 0.1, {'FontSize', 14});
% xtickangle(25);

%% Plot 2 - Correlation between channels (all features)
feature_index=find(contains(table2array(accum_mean_corr(:,1)), 'x'));
accum_table=accum_mean_corr(feature_index,:);
unique_type = unique(accum_table(:, 1), 'stable');
figure('pos',[250 250 1200 800]);

bar_data = [];
for i = 1:size(unique_type, 1)
   filter = accum_table(accum_table.Channel == table2array(unique_type(i, :)), :);

   bar_data = [bar_data; table2array(filter(:, 4))'];
end

bar_data = bar_data;
standard_errors = standard_errors;

%standard_error_reshape = reshape(standard_errors, 3, 6)';
barplot = barwitherr(standard_errors(feature_index)', bar_data);

% barplot = bar(bar_data);
% l = cell(1,3);
% l{1} = 'EEG';
% l{2} = 'EEG + EOG';
% l{3} = 'EEG + EOG + EMG';
% legend(barplot, l);

x_label=[ ...
 "EEG x LEOG", ...
 "EEG x EMG", ...
 "LEOG x EMG", ...
 "LEOG x REOG", ...
 "LEOG+REOG x LEOG-REOG", ...
 "LEOG+REOG x LEOG*REOG", ...
 "LEOG-REOG x LEOG*REOG"];
text(1:length(bar_data),bar_data,num2str(bar_data),'vert','bottom','horiz','center', 'FontSize', 14); 
title(strcat('Correlation of determination (R^2) between channels - ALL features (n=', num2str(size(all_dataset, 1)), ')'), 'FontSize', 16, 'Interpreter', 'tex');
grid on;
xlabel('Channel', 'FontSize', 15);
ylabel('Correlation of determination (R^2)', 'FontSize', 15, 'Interpreter', 'tex');
set(gca, 'xticklabel', x_label);
xtickangle(25);

%% Plot 3 - Correlation between channels (single features)
feature_index=find(contains(table2array(accum_mean_corr(:,1)), '#'));
accum_table=accum_mean_corr(feature_index,:);
unique_type = unique(accum_table(:, 1), 'stable');
figure('pos',[250 250 1200 800]);

bar_data = [];
for i = 1:size(unique_type, 1)
   filter = accum_table(accum_table.Channel == table2array(unique_type(i, :)), :);

   bar_data = [bar_data; table2array(filter(:, 4))'];
end

bar_data = bar_data;
standard_errors = standard_errors;

%standard_error_reshape = reshape(standard_errors, 3, 6)';
barplot = barwitherr(standard_errors(feature_index)', bar_data);

% barplot = bar(bar_data);
% l = cell(1,3);
% l{1} = 'EEG';
% l{2} = 'EEG + EOG';
% l{3} = 'EEG + EOG + EMG';
% legend(barplot, l);

x_label=[ ...
 "EEG x LEOG", ...
 "EEG x EMG", ...
 "LEOG x EMG", ...
 "LEOG x REOG", ...
 "LEOG+REOG x LEOG-REOG", ...
 "LEOG+REOG x LEOG*REOG", ...
 "LEOG-REOG x LEOG*REOG"];
text(1:length(bar_data),bar_data,num2str(bar_data),'vert','bottom','horiz','center', 'FontSize', 14); 
title(strcat('Mean Correlation of determination (R^2) between channels - single feature (n=', num2str(size(all_dataset, 1)), ')'), 'FontSize', 16, 'Interpreter', 'tex');
grid on;
xlabel('Channels', 'FontSize', 15);
ylabel('Mean Correlation of determination (R^2)', 'FontSize', 15, 'Interpreter', 'tex');
set(gca, 'xticklabel', x_label);
xtickangle(25);
barColorMap = jet(size(bar_data, 1));


%% Plot 4 - EOG correlation matrix

figure;
corr=eog_grand_corr_matrix./size(all_dataset, 1);

max_elem=size(corr, 1);
x = repmat(1:max_elem,max_elem,1);
y = x';
imagesc(corr);

t = cell(max_elem, max_elem);
for i=1:max_elem
    for j=1:max_elem
        t(i, j) = cellstr(sprintf('%0.3f', corr(i,j)));
    end
end

title(strcat('Correlation matrix between different EOG derivations (n=', num2str(size(all_dataset, 1)), ')'), 'FontSize', 20, 'Interpreter', 'tex');

xlabel('Channel derivations', 'FontSize', 16);
ylabel('Channel derivations', 'FontSize', 16);
x_label = ["LEOG", "REOG", "LEOG-REOG", "LEOG*REOG"];
set(gca, 'xticklabel', x_label, 'FontSize', 13);
set(gca, 'XTick', 1:1:5);
set(gca, 'yticklabel', x_label, 'FontSize', 13);
set(gca, 'YTick', 1:1:5);
text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 14, ...
    'FontWeight', 'bold');

caxis([0 1])

c = hot;
c = flipud(c);
colormap(c(5));
colorbar;

%% USING CONVOLUTION
%%
% mC = ceil(num_to_keep/2); % Distance to center of submatrix
% M = ones(num_to_keep);
% Aconv = conv2(CM,M); % Do the convolution.
% [~,minColIdx] = min(min(Aconv(1+mC:end-mC,1+mC:end-mC))); % Find column center with smallest sum
% [~,minRowIdx] = min(min(Aconv(1+mC:end-mC,minColIdx+mC),[],2)); % Find row center with smlest sum
% minRowIdx = minRowIdx+mC-1; % Convoluted matrix is larger than A
% minColIdx = minColIdx+mC-1; % Convoluted matrix is larger than A
% range = -mC+1:mC-1;
% B = CM(minRowIdx+range, minColIdx+range);
% 
% B;

%%
function  [ values ] = find_top_X_features(MAT, SUBSET)

    c=corr(MAT);
    N=size(c,1);

    c_r2=c.^2;
    ccols=c(:);
    ccols2=c_r2(:);
    [sc,idx]=sort(abs(ccols));
    m=[floor(idx./N)+1 mod(idx,N)];

    mod_0=find(m(:,2)==0);
    m(mod_0, 2) = N;

    m=m(1:2:end, :);

    % for i = 1:size(m,1)
    unique_features = [];
    total = 0;
    i = 1;
    while length(unique(unique_features)) < SUBSET
        unique_features = [unique_features, m(i,1)];
        unique_features = [unique_features, m(i,2)];

        total = total + abs(c(m(i,1), m(i,2)));
        i = i+1;
    end

    values=unique_features(1:SUBSET);

end


function [ values ] = the_optimal_method( CM , num_to_keep)
    %the_iterative_method Takes correlation matrix CM and number to keep, returns list of people who should be kicked out 

    N = size(CM,1);  

    clear model;  
    names = strseq('x',[1:N]);  
    model.varnames = names;  
    model.Q = sparse(CM); % Gurobi needs a sparse matrix as input  
    model.A = sparse(ones(1,N));  
    model.obj = zeros(1,N);  
    model.rhs = num_to_keep;  
    model.sense = '=';  
    model.vtype = 'B';

    gurobi_write(model, 'qp.mps');

    results = gurobi(model);
    

    values = results.x;

end
