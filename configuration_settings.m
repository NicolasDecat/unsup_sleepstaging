% This is a configuration setting files that shared among all matlab files.
% This file is to facilitate all the environment-specific variables that
% could be used so that all .m files have no hardcoded path.

%% Pre HCTSA
BLOCKEDFLOAD_DIR='/Volumes/Spaceship/SleepResearch/unsup_sleepstaging/library/blockedfloader';
DATA_DIR = '/Volumes/Seagate Expansion Drive/EDF_data/CCSHS/polysomnography/edfs';
%DATA_DIR = '/Volumes/Spaceship/ccshs_edf';
EDF_FILE = 'ccshs-trec-1800076.edf';

%% Shared (Common)
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';

%% Post HCTSA
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
HCTSA_DATA_DIR='/Volumes/Spaceship/Voss_Lucid/ME_N3_bipolar/';
%HCTSA_DATA_DIR='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800001_HCTSA/1800001_EMG_5chan';
%HCTSA_DATA_DIR='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800748_HCTSA';

% Changes these per run
WHICH_DATA = 439; % Which dataset
DATASET_NAME = '1800439';

OUTPUT_STATS_FILENAME = strcat('SUP_UNSUP_HCTSA_200_DS',num2str(WHICH_DATA),'.mat');
NUM_CHANNELS = 7;
NUM_CHANNELS_TO_RUN = [1];
%CONFIGURATIONS_TO_RUN = ["BALANCED_LABELED", "UNBALANCED_LABELED_A" , "UNBALANCED_LABELED_B"];
CONFIGURATIONS_TO_RUN = ["BALANCED_LABELED"];
EXPS_TO_RUN=[5:5];
SINGLE_FEATURE_RESULT=strcat('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_', DATASET_NAME, '_HCTSA/SINGLE_FEATURES_DS', num2str(WHICH_DATA), '.mat');
SELECT_TOP_200_FEATURES=[96,148,67,176,148,172,170,74,48,197,196,101,123,131,121,132,1,147,101,133,33,125,184,124,62,141,2,150,96,181,48,196,77,150,26,125,65,133,125,134,77,133,51,197,1,56,73,84,9,125,121,143,114,131,38,158,47,165,198,198,168,12,22,147,7,147,113,168,163,150,47,178,31,197,78,142,87,165,53,197,95,143,85,198,96,131,51,152,154,87,294,346,265,374,346,370,368,272,246,395,394,299,321,329,319,330,199,345,299,331,231,323,382,322,260,339,200,348,294,379,246,394,275,348,224,323,263,331,323,332,275,331,249,395,199,254,271,282,207,323,319,341,312,329,236,356,245,363,396,396,366,210,220,345,205,345,311,366,361,348,245,376,229,395,276,340,285,363,251,395,293,341,283,396,294,329,249,350,352,285,339,388,337,343,310,338,278,395,348,393,492,544,463,572,544,568,566,470,444,593];


%% selectdata and crossval
ANSWER_FILE=strcat('/Volumes/Spaceship/ccshs_datasets/ccshs_', DATASET_NAME, '_annot.mat');
CONFUSION_MATRIX_FILE=strcat('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_', DATASET_NAME, '_HCTSA');
% CM_SAVE_DIR='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800001_HCTSA';
% HCTSA_FILE='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800001_HCTSA/HCTSA_N.mat';
CM_SAVE_DIR=strcat('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_', DATASET_NAME, '_HCTSA');
HCTSA_FILE=strcat('/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_', DATASET_NAME, '_HCTSA/HCTSA_N.mat');
OPS_FILE='reduced_ops.txt';
TRAINING_PERCENTAGE = 0.7;
NUM_CLUSTERS = 5;
CROSSVAL_ITERATION=5;
%CUSTOM_CHANNELS = [1374*0+1:1374*1];
PLOT_CONFUSION_MATRIX=true;
PLOT_ACCURACY_REPORT=false;

%% Experiment configurations to run
THREAD_POOL=4;
EVAL_THRESHOLD=0.05;
TOP_FEATURE_COUNT=200;
FSLIB_TOOLBOX_DIR='/Users/Zhao/Documents/MATLAB/Add-Ons/Toolboxes/Feature Selection Library/code/FSLib_v5.2_2017';
% This will output the data and images for each iteration. This requires the image dir below.
DEBUG_CROSSVALIDATION=false; 
DEBUG_CROSSVALIDATION_IMAGEDIR='/Volumes/Spaceship/ccshs_1800007_3EXG';
