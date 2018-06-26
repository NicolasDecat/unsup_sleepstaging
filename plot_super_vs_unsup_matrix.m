configuration_settings;

homedir = pwd;
%% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)

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


t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/BALANCED_LABELED_EEG_1_AVG.mat');
eeg_balanced_labeled_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/BALANCED_LABELED_EOG_1_AVG.mat');
eog_balanced_labeled_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/BALANCED_LABELED_EMG_1_AVG.mat');
emg_balanced_labeled_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/UNBALANCED_LABELED_A_EEG_2_AVG.mat');
eeg_unbalanced_labeled_A_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/UNBALANCED_LABELED_A_EOG_2_AVG.mat');
eog_unbalanced_labeled_A_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/UNBALANCED_LABELED_A_EMG_2_AVG.mat');
emg_unbalanced_labeled_A_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/UNBALANCED_LABELED_B_EEG_3_AVG.mat');
eeg_unbalanced_labeled_B_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/UNBALANCED_LABELED_B_EOG_3_AVG.mat');
eog_unbalanced_labeled_B_stats=t.save_stats;

t=load('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/UNBALANCED_LABELED_B_EMG_3_AVG.mat');
emg_unbalanced_labeled_B_stats=t.save_stats;

%% Define the matrix
dm = [ ...
    table2array(eeg_balanced_labeled_stats(eeg_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", 6)), ...
    table2array(eeg_unbalanced_labeled_A_stats(eeg_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", 6)), ...
    table2array(eeg_unbalanced_labeled_B_stats(eeg_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", 6)), ...
    table2array(eeg_balanced_labeled_stats(eeg_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", 6)), ...
    table2array(eeg_unbalanced_labeled_A_stats(eeg_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", 6)), ...
    table2array(eeg_unbalanced_labeled_B_stats(eeg_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", 6)); ...
    table2array(eog_balanced_labeled_stats(eog_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", 6)), ...
    table2array(eog_unbalanced_labeled_A_stats(eog_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", 6)), ...
    table2array(eog_unbalanced_labeled_B_stats(eog_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", 6)), ...
    table2array(eog_balanced_labeled_stats(eog_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", 6)), ...
    table2array(eog_unbalanced_labeled_A_stats(eog_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", 6)), ...
    table2array(eog_unbalanced_labeled_B_stats(eog_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", 6)); ...
    table2array(emg_balanced_labeled_stats(emg_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", 6)), ...
    table2array(emg_unbalanced_labeled_A_stats(emg_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", 6)), ...
    table2array(emg_unbalanced_labeled_B_stats(emg_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", 6)), ...
    table2array(emg_balanced_labeled_stats(emg_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", 6)), ...    
    table2array(emg_unbalanced_labeled_A_stats(emg_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", 6)), ...
    table2array(emg_unbalanced_labeled_B_stats(emg_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", 6)) ...
    ];

x_label = ["Supervised Balanced Labeled", "Supervised Unbalanced Labeled A", "Supervised Unbalanced Labeled B", ...
    "Unsupervised Balanced Labeled",  "Unsupervised Unbalanced Labeled A", "Unsupervised Unbalanced Labeled B"];

no_features=size(unique(eeg_balanced_labeled_stats(:,2)),1);
assert(sum(table2array(unique(eeg_balanced_labeled_stats(:,2))) == table2array(unique(eog_balanced_labeled_stats(:,2))))'==no_features);
assert(sum(table2array(unique(eeg_balanced_labeled_stats(:,2))) == table2array(unique(emg_balanced_labeled_stats(:,2))))'==no_features);
assert(sum(table2array(unique(eeg_unbalanced_labeled_A_stats(:,2))) == table2array(unique(eog_unbalanced_labeled_A_stats(:,2))))'==no_features);
assert(sum(table2array(unique(eeg_unbalanced_labeled_A_stats(:,2))) == table2array(unique(emg_unbalanced_labeled_A_stats(:,2))))'==no_features);
assert(sum(table2array(unique(eeg_unbalanced_labeled_B_stats(:,2))) == table2array(unique(eog_unbalanced_labeled_B_stats(:,2))))'==no_features);
assert(sum(table2array(unique(eeg_unbalanced_labeled_B_stats(:,2))) == table2array(unique(emg_unbalanced_labeled_B_stats(:,2))))'==no_features);

y_label = [strcat('EEG-', table2array(unique(eeg_balanced_labeled_stats(:,2)))); ...
    strcat('EOG-', table2array(unique(eog_balanced_labeled_stats(:,2)))); ...
    strcat('EMG-', table2array(unique(emg_balanced_labeled_stats(:,2))))];


%% Draw figure
figure;
title('Features by configurations');

%subplot(2, 6, [1 2 3 4 5 6]);
imagesc(dm);
colorbar('Location', 'westoutside');

customColorMap = flipud(BF_getcmap('redyellowblue',6,0));
colormap(customColorMap);

ax=gca;
ax.XAxis.TickValues = [1 2 3 4 5 6];
xtickangle(0);
ax.XAxis.TickLabels = x_label;
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%% Draw another figure
figure;
colorbar('Location', 'westoutside');
colormap(customColorMap);

%
subplot(2, 6, [1 7]);

e1=eeg_balanced_labeled_stats(eeg_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", :);
e2=eog_balanced_labeled_stats(eog_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", :);
e3=emg_balanced_labeled_stats(emg_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", :);
e1=sortrows(e1, 'TestingAccuracy', 'descend');
e2=sortrows(e2, 'TestingAccuracy', 'descend');
e3=sortrows(e3, 'TestingAccuracy', 'descend');
y_label = [strcat('EEG-', table2array(unique(e1(:,2)))); ...
    strcat('EOG-', table2array(e2(:,2))); ...
    strcat('EMG-', table2array(e3(:,2)))];

imagesc([double(table2array(e1(:, 6))); double(table2array(e2(:, 6))); double(table2array(e3(:, 6)))]);
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(1);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [2 8]);

e1=eeg_unbalanced_labeled_A_stats(eeg_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", :);
e2=eog_unbalanced_labeled_A_stats(eog_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", :);
e3=emg_unbalanced_labeled_A_stats(emg_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", :);
e1=sortrows(e1, 'TestingAccuracy', 'descend');
e2=sortrows(e2, 'TestingAccuracy', 'descend');
e3=sortrows(e3, 'TestingAccuracy', 'descend');
y_label = [strcat('EEG-', table2array(unique(e1(:,2)))); ...
    strcat('EOG-', table2array(e2(:,2))); ...
    strcat('EMG-', table2array(e3(:,2)))];

imagesc([double(table2array(e1(:, 6))); double(table2array(e2(:, 6))); double(table2array(e3(:, 6)))]);
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(2);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [3 9]);

e1=eeg_unbalanced_labeled_B_stats(eeg_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", :);
e2=eog_unbalanced_labeled_B_stats(eog_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", :);
e3=emg_unbalanced_labeled_B_stats(emg_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", :);
e1=sortrows(e1, 'TestingAccuracy', 'descend');
e2=sortrows(e2, 'TestingAccuracy', 'descend');
e3=sortrows(e3, 'TestingAccuracy', 'descend');
y_label = [strcat('EEG-', table2array(unique(e1(:,2)))); ...
    strcat('EOG-', table2array(e2(:,2))); ...
    strcat('EMG-', table2array(e3(:,2)))];

imagesc([double(table2array(e1(:, 6))); double(table2array(e2(:, 6))); double(table2array(e3(:, 6)))]);
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(3);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [4 10]);

e1=eeg_balanced_labeled_stats(eeg_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", :);
e2=eog_balanced_labeled_stats(eog_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", :);
e3=emg_balanced_labeled_stats(emg_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", :);
e1=sortrows(e1, 'TestingAccuracy', 'descend');
e2=sortrows(e2, 'TestingAccuracy', 'descend');
e3=sortrows(e3, 'TestingAccuracy', 'descend');
y_label = [strcat('EEG-', table2array(unique(e1(:,2)))); ...
    strcat('EOG-', table2array(e2(:,2))); ...
    strcat('EMG-', table2array(e3(:,2)))];

imagesc([double(table2array(e1(:, 6))); double(table2array(e2(:, 6))); double(table2array(e3(:, 6)))]);
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(4);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [5 11]);

e1=eeg_unbalanced_labeled_A_stats(eeg_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", :);
e2=eog_unbalanced_labeled_A_stats(eog_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", :);
e3=emg_unbalanced_labeled_A_stats(emg_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", :);
e1=sortrows(e1, 'TestingAccuracy', 'descend');
e2=sortrows(e2, 'TestingAccuracy', 'descend');
e3=sortrows(e3, 'TestingAccuracy', 'descend');
y_label = [strcat('EEG-', table2array(unique(e1(:,2)))); ...
    strcat('EOG-', table2array(e2(:,2))); ...
    strcat('EMG-', table2array(e3(:,2)))];

imagesc([double(table2array(e1(:, 6))); double(table2array(e2(:, 6))); double(table2array(e3(:, 6)))]);
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(5);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [6 12]);

e1=eeg_unbalanced_labeled_B_stats(eeg_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", :);
e2=eog_unbalanced_labeled_B_stats(eog_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", :);
e3=emg_unbalanced_labeled_B_stats(emg_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", :);
e1=sortrows(e1, 'TestingAccuracy', 'descend');
e2=sortrows(e2, 'TestingAccuracy', 'descend');
e3=sortrows(e3, 'TestingAccuracy', 'descend');
y_label = [strcat('EEG-', table2array(unique(e1(:,2)))); ...
    strcat('EOG-', table2array(e2(:,2))); ...
    strcat('EMG-', table2array(e3(:,2)))];

imagesc([double(table2array(e1(:, 6))); double(table2array(e2(:, 6))); double(table2array(e3(:, 6)))]);
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(6);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%% Draw another figure figure

%% Draw another figure
figure1=figure;
colorbar('Location', 'westoutside');
customColorMap = flipud(BF_getcmap('redyellowblue',3,0));
colormap(customColorMap);

%
subplot(2, 6, [1 7]);

e1=eeg_balanced_labeled_stats(eeg_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", :);
e1.FeatureId=strcat('EEG-', e1.FeatureId);
e2=eog_balanced_labeled_stats(eog_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", :);
e2.FeatureId=strcat('EOG-', e2.FeatureId);
e3=emg_balanced_labeled_stats(emg_balanced_labeled_stats.Type=="Supervised_Balanced_Labeled", :);
e3.FeatureId=strcat('EMG-', e3.FeatureId);

e=vertcat(e1,e2,e3);
e=sortrows(e, 'TestingAccuracy', 'descend');
save_stats=e;
save(strcat(CM_SAVE_DIR, filesep, 'SUPERVISED_BALANCED_LABELED_FEATURES_ORDERED.mat'), 'save_stats');

idx=find(e.NumberOfChannels=="EEG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*3);
idx=find(e.NumberOfChannels=="EOG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*2);
idx=find(e.NumberOfChannels=="EMG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*1);

y_label = table2array(e(:,2));

imagesc(double(table2array(e(:, 4))));
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(1);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [2 8]);

e1=eeg_unbalanced_labeled_A_stats(eeg_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", :);
e1.FeatureId=strcat('EEG-', e1.FeatureId);
e2=eog_unbalanced_labeled_A_stats(eog_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", :);
e2.FeatureId=strcat('EOG-', e2.FeatureId);
e3=emg_unbalanced_labeled_A_stats(emg_unbalanced_labeled_A_stats.Type=="Supervised_Unbalanced_Labeled_A", :);
e3.FeatureId=strcat('EMG-', e3.FeatureId);

e=vertcat(e1,e2,e3);
e=sortrows(e, 'TestingAccuracy', 'descend');
save_stats=e;
save(strcat(CM_SAVE_DIR, filesep, 'SUPERVISED_UNBALANCED_LABELED_A_FEATURES_ORDERED.mat'), 'save_stats');

idx=find(e.NumberOfChannels=="EEG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*3);
idx=find(e.NumberOfChannels=="EOG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*2);
idx=find(e.NumberOfChannels=="EMG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*1);

y_label = table2array(e(:,2));
imagesc(double(table2array(e(:, 4))));
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(2);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [3 9]);

e1=eeg_unbalanced_labeled_B_stats(eeg_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", :);
e1.FeatureId=strcat('EEG-', e1.FeatureId);
e2=eog_unbalanced_labeled_B_stats(eog_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", :);
e2.FeatureId=strcat('EOG-', e2.FeatureId);
e3=emg_unbalanced_labeled_B_stats(emg_unbalanced_labeled_B_stats.Type=="Supervised_Unbalanced_Labeled_B", :);
e3.FeatureId=strcat('EMG-', e3.FeatureId);

e=vertcat(e1,e2,e3);
e=sortrows(e, 'TestingAccuracy', 'descend');
save_stats=e;
save(strcat(CM_SAVE_DIR, filesep, 'SUPERVISED_UNBALANCED_LABELED_B_FEATURES_ORDERED.mat'), 'save_stats');

idx=find(e.NumberOfChannels=="EEG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*3);
idx=find(e.NumberOfChannels=="EOG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*2);
idx=find(e.NumberOfChannels=="EMG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*1);

imagesc(double(table2array(e(:, 4))));
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(3);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [4 10]);

e1=eeg_balanced_labeled_stats(eeg_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", :);
e1.FeatureId=strcat('EEG-', e1.FeatureId);
e2=eog_balanced_labeled_stats(eog_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", :);
e2.FeatureId=strcat('EOG-', e2.FeatureId);
e3=emg_balanced_labeled_stats(emg_balanced_labeled_stats.Type=="Unsupervised_Balanced_Labeled", :);
e3.FeatureId=strcat('EMG-', e3.FeatureId);

e=vertcat(e1,e2,e3);
e=sortrows(e, 'TestingAccuracy', 'descend');
save_stats=e;
save(strcat(CM_SAVE_DIR, filesep, 'UNSUPERVISED_BALANCED_LABELED_FEATURES_ORDERED.mat'), 'save_stats');

idx=find(e.NumberOfChannels=="EEG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*3);
idx=find(e.NumberOfChannels=="EOG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*2);
idx=find(e.NumberOfChannels=="EMG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*1);

y_label = table2array(e(:,2));
imagesc(double(table2array(e(:, 4))));
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(4);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [5 11]);

e1=eeg_unbalanced_labeled_A_stats(eeg_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", :);
e1.FeatureId=strcat('EEG-', e1.FeatureId);
e2=eog_unbalanced_labeled_A_stats(eog_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", :);
e2.FeatureId=strcat('EOG-', e2.FeatureId);
e3=emg_unbalanced_labeled_A_stats(emg_unbalanced_labeled_A_stats.Type=="Unsupervised_Unbalanced_Labeled_A", :);
e3.FeatureId=strcat('EMG-', e3.FeatureId);

e=vertcat(e1,e2,e3);
e=sortrows(e, 'TestingAccuracy', 'descend');
save_stats=e;
save(strcat(CM_SAVE_DIR, filesep, 'UNSUPERVISED_UNBALANCED_LABELED_A_FEATURES_ORDERED.mat'), 'save_stats');

idx=find(e.NumberOfChannels=="EEG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*3);
idx=find(e.NumberOfChannels=="EOG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*2);
idx=find(e.NumberOfChannels=="EMG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*1);

y_label = table2array(e(:,2));
imagesc(double(table2array(e(:, 4))));
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(5);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

%

subplot(2, 6, [6 12]);

e1=eeg_unbalanced_labeled_B_stats(eeg_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", :);
e1.FeatureId=strcat('EEG-', e1.FeatureId);
e2=eog_unbalanced_labeled_B_stats(eog_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", :);
e2.FeatureId=strcat('EOG-', e2.FeatureId);
e3=emg_unbalanced_labeled_B_stats(emg_unbalanced_labeled_B_stats.Type=="Unsupervised_Unbalanced_Labeled_B", :);
e3.FeatureId=strcat('EMG-', e3.FeatureId);

e=vertcat(e1,e2,e3);
e=sortrows(e, 'TestingAccuracy', 'descend');
save_stats=e;
save(strcat(CM_SAVE_DIR, filesep, 'UNSUPERVISED_UNBALANCED_LABELED_B_FEATURES_ORDERED.mat'), 'save_stats');

idx=find(e.NumberOfChannels=="EEG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*3);
idx=find(e.NumberOfChannels=="EOG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*2);
idx=find(e.NumberOfChannels=="EMG");
e(idx, 4) = array2table(ones(size(idx,1), 1)*1);

y_label = table2array(e(:,2));
imagesc(double(table2array(e(:, 4))));
ax=gca;
xtickangle(0);
ax.XAxis.TickValues = [1];
ax.XAxis.TickLabels = x_label(6);
ax.YAxis.TickValues = [1:size(y_label, 1)];
ax.YAxis.TickLabels = y_label';

colorbar('Ticks',[1 2 3],'TickLabels',{'EMG','EOG','EEG'},'Position',[0.92265625 0.11031175059952 0.01484375 0.815347721822542]);


%% Display the top selected features in CSV.
% EEG_TOP = 100;
% EOG_TOP = 80;
% EMG_TOP = 20;
files = ["SUPERVISED_BALANCED_LABELED_FEATURES_ORDERED.mat", ...
    "SUPERVISED_UNBALANCED_LABELED_A_FEATURES_ORDERED.mat", ...
    "SUPERVISED_UNBALANCED_LABELED_B_FEATURES_ORDERED.mat", ...
    "UNSUPERVISED_BALANCED_LABELED_FEATURES_ORDERED.mat", ...
    "UNSUPERVISED_UNBALANCED_LABELED_A_FEATURES_ORDERED.mat", ...
    "UNSUPERVISED_UNBALANCED_LABELED_B_FEATURES_ORDERED.mat"];

for i=1:length(files)
   file = files(i);
   
   ftable = load(strcat(CM_SAVE_DIR, filesep, file));
   ftable = ftable.save_stats;
   
   eeg_features = ftable(ftable.NumberOfChannels=="EEG", [2, 3, 6]);
   eog_features = ftable(ftable.NumberOfChannels=="EOG", [2, 3, 6]);
   emg_features = ftable(ftable.NumberOfChannels=="EMG", [2, 3, 6]);
%    eeg_eog_features = ftable(ftable.NumberOfChannels=="EEG" | ftable.NumberOfChannels=="EOG", [2, 3, 6]);
   all_features = ftable(:, [2, 3, 6]);
   
   %featureMat = strings(EEG_TOP+EOG_TOP+EMG_TOP, 2);
%    featureMat = [table2array(eeg_features(1:EEG_TOP, [1, 2])); table2array(eog_features(1:EOG_TOP, [1, 2])); ...
%        table2array(emg_features(1:EMG_TOP, [1, 2]))];
   featureMat = [table2array(all_features(1:198, [1, 2]))];
   
   writetable(array2table(featureMat), strcat(CM_SAVE_DIR, filesep, file, "_TOP_200_EEG_EOG_EMG_FEATURES.csv"));
   
end
%%
set(gcf, 'units', 'normalized', 'position', [0.05 0.15 0.9 0.8])





