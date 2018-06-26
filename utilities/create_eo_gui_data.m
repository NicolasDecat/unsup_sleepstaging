clc;

%% 
% This is a script that prepare data for EO tool to analyse eye movement.

%% CONFIGURATION
EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/KJ_N2_data.mat';
%EDF_FILE='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging_master/Voss_Pilot/SK_N1_F1_data.mat';
%EDF_FILE='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging_master/Voss_Pilot/SK_N1_F1_data.mat';

OUTPUT_FILE='KJ_N2_EOG.mat';

% Global

load(EDF_FILE);

%% Extract channel names, sampling frequency from header struct
channel = cellstr(header.channelname);  % Channel names + converted into cell array
fs = header.samplerate(1);      % Assume sampling rate are the same for all channels
recordtime = header.records;    % aaTotal recorded time (seconds)
nchannels = header.channels;    % Number of channels

channel_info=[(1:length(header.channelname))',string(header.channelname)]

% Referencing
% data(gci(channel_info, "F3"),:)=data(gci(channel_info, "F3"),:) - data(gci(channel_info, "Cz"), :);
EOG_horizontal_channel=gci(channel_info, "EOGright");
EOG_vertical_channel=gci(channel_info, "EOGup");

%%
samplingrate=fs; %sampling rate in Hertz
mspersample=1000/samplingrate; %milliseconds per sample 
eodata=data([EOG_horizontal_channel EOG_vertical_channel],:);
eodata=eodata';
mytime=[1:1:size(eodata,1)]'; %create column vector of length same as your data 
mytime=mytime*mspersample; %adjust time vector to sample rate
data= [mytime eodata]; %attach the time column to your x and y data already stored in 'data'

save(strcat("/Users/Zhao/Downloads", filesep, OUTPUT_FILE), 'data');


function index = gci(channel_info, chn_name) 
    index = find(strcmp(strtrim(channel_info(:,2)), chn_name));
end