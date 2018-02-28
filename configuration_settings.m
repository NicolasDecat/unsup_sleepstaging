% This is a configuration setting files that shared among all matlab files.
% This file is to facilitate all the environment-specific variables that
% could be used so that all .m files have no hardcoded path.

% Pre HCTSA
BLOCKEDFLOAD_DIR='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/library/blockedfloader';
DATA_DIR = '/Volumes/Untitled/edfs';
EDF_FILE = 'ccshs-trec-1800001.edf'
WHICH_DATA = 1;

% Shared (Common)
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
NUM_CHANNELS = 3;

% Post HCTSA
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa';
HCTSA_DATA_DIR='./180001_HCTSA';

% selectdata and crossval
ANSWER_FILE='ccshs_1800001_annot.mat';
CONFUSION_MATRIX_FILE='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/180001_HCTSA';
CM_SAVE_DIR='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/180001_HCTSA';
HCTSA_FILE='180001_HCTSA/HCTSA_N.mat';
OPS_FILE='reduced_ops.txt';
TRAINING_PERCENTAGE = 0.7;
NUM_CHANNELS_USED_FOR_CROSSVAL = 1;
FSLIB_TOOLBOX_DIR='/Users/Zhao/Documents/MATLAB/Add-Ons/Toolboxes/Feature Selection Library/code/FSLib_v5.2_2017';
