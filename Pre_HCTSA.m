%% Reading EDF format into MAT using EEGLAB
% EEGLAB: MATLAB toolbox https://sccn.ucsd.edu/eeglab/
basedir = pwd;
edffile = ('learn-nsrr01.edf'); % Change filename for new data set

%% Add eeglab to path to use readedf function
addpath(genpath('C:\Users\Piengkwan\Documents\MATLAB\eeglab'))  % Change directory accordingly
% Read .edf data into workspace and save as .mat
% Extract channels information, sampling frequency

cd('.\310817\learn\polysomnography\edfs')

% Possible to read without initiating EEGLAB?
% [data,header]=readedf(edffile); 
cd(basedir)
%% Load edf file
eeglab
EEG = pop_biosig(['C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\310817\learn\polysomnography\edfs\',edffile]);

%% Load data manually (Only for 'Learn' dataset because the naming
% convention is different - cannot use readedf() directly)
data = EEG.data;
header.channels = EEG.nbchan;
header.records = EEG.pnts;
header.samplerate = EEG.srate;
for n = 1:EEG.nbchan
    header.channelname(n,1) = cellstr(EEG.chanlocs(n).labels);
end

% Save data file
dataname = sprintf([edffile(1:length(edffile)-4),'_','data']);
save(dataname,'data','header')
% save in base directory, for further use or modification of time segment
% length
%% After converting from .edf to .mat
cd('.\010917_Learn01')
addpath(genpath('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging'))
load('learn-nsrr01_data.mat')
% Extract channel names, sampling frequency from header struct
channel = cellstr(header.channelname);  % Channel names + converted into cell array
fs = header.samplerate(1);      % Assume sampling rate are the same for all channels
recordtime = header.records;    % Total recorded time (seconds)
nchannels = header.channels;    % Number of channels

% Choose only EEG, EMG, EOG and ECG channels
indexE_G = [3,5,6,7,8];
channelE_G = channel(indexE_G);
dataE_G = data(indexE_G,:);
nchanE_G = length(indexE_G);
% Read data and segment into specified length 
% 30-second epoch - For comparison with human performance)
% 5-second epoch - For detection substages
interval = 30;
timeSeriesData = read_edf_segment(dataE_G,fs,interval,nchanE_G);% Change function name


% Number of time segment/epoch
n_seg = (recordtime/fs)/interval;
% n_seg = size(timeSeriesData,2);
n_ts = size(timeSeriesData,1);  % Number of row of timeseriesData

% Time segment name
for i=1:n_seg
    name = sprintf('timeseg_%d',i);
    timelabel{i}=name;
end

[labels,keywords] = lagelgen(n_ts,2,channelE_G,timelabel);

% Save in HCTSA input format
hctsafile = sprintf(['TS_',edffile(1:length(edffile)-4),'_',num2str(interval),'sec']);
save(hctsafile,'timeSeriesData','labels','keywords')

