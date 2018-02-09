%% Reading .edf data and convert into .mat
% The EEG data was in EDF format, this code is the first step before
% Pre_HCTSA.m
% The readedf function is in EEGLAB toolbox (https://sccn.ucsd.edu/eeglab/)

basedir = pwd;
% cd (location of edf files)
edffile = ('ME_N3.edf'); % Change filename for new data set

% Add eeglab to path to use readedf function
addpath(genpath('C:\Users\Piengkwan\Documents\MATLAB\eeglab'))
% Read data into mattrix
[data,header]=readedf(edffile); 

% Name mat file according to edf filename
dataname = sprintf([edffile(1:length(edffile)-4),'_','data']);
mkdir 041017_RawData    % Change name of folder to store data in 
% addpath(genpath('.\041017_RawData'))
cd ./041017_RawData
save(dataname,'data','header')
% save for further use or modification of time segment length
cd(basedir)
