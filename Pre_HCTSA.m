configuration_settings

%% Preprocessing of polysomnographic data and generate time-series matrix for HCTSA
homedir = pwd;
edffile = strcat(EDF_FOLDER, filesep, EDF_FILE); 
% load(edffile)
cd(homedir)
%% Load data using blockEdfLoad
% https://sleepdata.org/tools/dennisdean-block-edf-loader
addpath(genpath(BLOCKEDFLOAD_DIR));
[header, signalHeader, signalCell] = blockEdfLoad(edffile);

%% Extract channel names, sampling frequency for 1 EEG re-ref channel
channel = strcat(signalHeader(1).signal_labels,'-',signalHeader(4).signal_labels);
fs = [signalHeader(1).samples_in_record, signalHeader(4).samples_in_record]; % Sampling rate of the C3 channel
recordtime = header.num_data_records; % Total recording time (seconds)
nchannels = 1; % Only 1 EEG channel
% Check if the bipolar re-reference pair is valid
if fs(1)~=fs(2)
    valid = 0;
else
    valid = 1;
    fs = fs(1);
end
%% (Additional part)Pre-processing
% Bipolar re-reference C3-A2 (Channel derivation from CCSHS protocol)
data = transpose(signalCell{1}-signalCell{4});
%% Read data and segment into specified length 
% 25-second epoch - For comparison with human performance)
% 5-second epoch - For detection substages
interval = 30;
timeSeriesData = read_edf_segment(data,fs,interval,nchannels);% Change function name

% Number of time segment/epoch
n_seg = recordtime; % Record time (second)
[n_ts,~] = size(timeSeriesData);
%% Time segment name
for i=1:n_ts
    name = sprintf('timeseg_%d',i)
    timelabel{i}=name;
end
%%
[labels,keywords] = labelgen(n_ts,1,timelabel);

%% Save in HCTSA input format
hctsafile = sprintf(['TS_',edffile(1:length(edffile)-4)]);
% mkdir 041017_KJ_N2
% cd ./041017_KJ_N2
save(hctsafile,'timeSeriesData','labels','keywords')

