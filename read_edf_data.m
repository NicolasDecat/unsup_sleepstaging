% Whole night recording, all channels(EEG,EMG,EOG,ECG)

%% Load data into workspace after edf is converted into .mat with 'data' and 'header' variables
load('KJ_N1_data.mat')

% Get channel name from header for labelling time series  
chanlabel = cellstr(header.channelname);

%% Segment into 25-sec interval (5000 samples)
samplingF = 200; % Hz
interval = 25; % (seconds) specified, use 25 or 5 for substaging 
channel = 26; % number of channels in this set of data

timeSeriesData = read_edf_segment(data,samplingF,interval,channel); % Segment using read_edf_segment function

%%  Generate labels and keywords
n_ts = length(timeSeriesData);
%%  Level 2 labels
% Number of time segments
nSeg = 
for n = 1:948 % change to variable 
    name = sprintf('timeseg_%d',n);
    level2{n} = name;
end
%% Label and keyowrds
[labels,keywords] = labelgen(n_ts,2,chanlabel,level2,level3); % Need to edit label function again

%% Save file + ready for HCTSA
save('TS_KJ_N1_25sec.mat','timeSeriesData','labels','keywords')