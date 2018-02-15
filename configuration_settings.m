% This is a configuration setting files that shared among all matlab files.
% This file is to facilitate all the environment-specific variables that
% could be used so that all .m files have no hardcoded path.

% Pre HCTSA
%BLOCKEDFLOAD_DIR='C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\blockEdfLoad'
BLOCKEDFLOAD_DIR='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/library/blockedfloader'
EDF_FOLDER = '/Volumes/Untitled/edfs'
EDF_FILE = 'ccshs-trec-1800032.edf'

% Post HCTSA
HCTSA_DATA_DIR='./180001_HCTSA'
HCTSA_DIR='/Users/Zhao/SleepPsychoPhysics/Source/hctsa'

% selectdata and crossval
ANSWER_FILE='ccshs_1800001_annot.mat'
CONFUSION_MATRIX_FILE='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/180001_HCTSA'
CM_SAVE_DIR='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging/180001_HCTSA'
HCTSA_FILE='180001_HCTSA/HCTSA_N.mat'
REDUCE_OPTS_FILE='reduced_ops.txt'

