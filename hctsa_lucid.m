%% Configuration

HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
TARGET_FOLDER='/Volumes/Spaceship/Voss_Lucid/KJ_N1/ALL_EEG';
%TARGET_FILE='HCTSA_N_1_EEG_2_substage_Full.mat';
%TARGET_FILE='HCTSA_N_1_EEG_2_substage_Sub.mat';
%TARGET_FILE='HCTSA_N_1_EEG_4_substage_Full.mat';
%TARGET_FILE='HCTSA_N_1_EEG_4_substage_Sub.mat';
%TARGET_FILE='HCTSA_N_3_EEG_2_substage_Full.mat';
%TARGET_FILE='HCTSA_N_3_EEG_2_substage_Sub.mat';
%TARGET_FILE='HCTSA_N_3_EEG_4_substage_Full.mat';
%TARGET_FILE='HCTSA_N_3_EEG_4_substage_Sub.mat';
%TARGET_FILE='HCTSA_N_3_EEG_3_substage_Full.mat';
%TARGET_FILE='HCTSA_N_3_EEG_3_substage_Sub.mat';
%TARGET_FILE='HCTSA_N_1_EEG_Frontal_2_substage_Full.mat';
TARGET_FILE='HCTSA_N_1_EEG_Frontal_2_substage_Sub.mat';

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

%%
cd(TARGET_FOLDER)

% Copy the current file to standard HCTSA.mat naming for LabelGroups
copyfile(TARGET_FILE, 'HCTSA.mat');

TS_LabelGroups([]);

movefile('HCTSA.mat', 'HCTSA_N.mat');
TS_cluster(distanceMetricRow, linkageMethodRow, distanceMetricCol, linkageMethodCol);

%TS_plot_DataMatrix('norm');

TS_plot_DataMatrix('cl', 'colorGroups', false);
TS_plot_DataMatrix('cl', 'colorGroups', true);

cd(homedir);
