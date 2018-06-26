%% Configuration

HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
TARGET_FOLDER='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA_Analysis';
TARGET_FILE='HCTSA_180001_N.mat';
ANSWER_FILE='/Volumes/Spaceship/ccshs_datasets/ccshs_1800001_annot.mat';

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

annotation = load(ANSWER_FILE);
TS_LabelGroups([]);

movefile('HCTSA.mat', 'HCTSA_N.mat');
TS_cluster(distanceMetricRow, linkageMethodRow, distanceMetricCol, linkageMethodCol);

%TS_plot_DataMatrix('norm');

TS_plot_DataMatrix('cl', 'colorGroups', false);
TS_plot_DataMatrix('cl', 'colorGroups', true);

cd(homedir);
