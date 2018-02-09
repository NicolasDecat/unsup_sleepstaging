%% Pre_HCTSA.m for labeled data
addpath(genpath('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging'))

% Load .mat file containing polysognograms data
basedir = pwd;
edffile = ('learn-nsrr02.edf'); % Change filename for new data set
cd('.\010917_Learn02')
load('learn-nsrr02_data.mat')
cd(basedir)

% Extract channel names, sampling frequency from header struct
channel = cellstr(header.channelname);  % Channel names + converted into cell array
fs = header.samplerate(1);      % Assume sampling rate are the same for all channels
recordtime = header.records;    % Total recorded time (seconds)
nchannels = header.channels;    % Number of channels

% Choose only wanteed channels (EEG, EOG, EMG)
indexE_G = [3,5,6,7,8]; % EEG: [3,8]| ExG:[3,5,6,7,8];
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
%%
[labels,keywords] = lagelgen(n_ts,2,channelE_G,timelabel);

%% Add scored sleep stages as keywords
% Load data
cd('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\310817\learn\polysomnography\annotations-events-profusion')
load('Learn02Annot.mat')
cd(basedir)
% sleepstage
%% Generate Keywords in TimeSeries
% Sleep stage for each channel
sleepStg = repmat(sleepstage,[nchanE_G,1]);
%%
stageName = {'W','N1','N2','N3','N4','REM'};
%%
stagelabel = cell(size(sleepStg));
for i=1:length(sleepStg)
    StgName = stageName(sleepStg(i)+1);
    stagelabel(i)= {StgName};
end
%% Append to existing keywords
for j=1:length(stagelabel)
    oldKeyword = cellstr(keywords(j));
    newKeyword = string(stagelabel(j));
    kkeywordss(j) = cellstr(strcat(oldKeyword,',',newKeyword));
 end
%% Save in HCTSA input format
keywords = kkeywordss;
cd('081017_Learn02_Labeled')
hctsafile = sprintf(['TS_',edffile(1:length(edffile)-4),'_',num2str(interval),'sec']);
save(hctsafile,'timeSeriesData','labels','keywords')


% EEG channel: [3,8]