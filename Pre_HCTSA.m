%% Preprocessing of EEG data and generate time-series matrix for HCTSA
homedir = pwd;
cd 041017_Rawdata
edffile = 'KJ_N2_data.mat'; 
load(edffile)
cd(homedir)
%% Extract channel names, sampling frequency from header struct
channel = cellstr(header.channelname);  % Channel names + converted into cell array
fs = header.samplerate(1);      % Assume sampling rate are the same for all channels
recordtime = header.records;    % Total recorded time (seconds)
nchannels = header.channels;    % Number of channels

%% (Additional part)
% For preprocessing and channel selection
%% Select 6 channels (EOG vertical, EOG horizontal, EMG, Frontal, Parietal, Occipital)
% EOG vertical = EOGup - EOGdown (?)
% figure;
% up = subplot(4,1,1);
% plot(data(1,1:120000))
% legend('EOGup')
% down = subplot(4,1,2);
% plot(data(2,1:120000))
% legend('EOGdown')
% plus = subplot(4,1,3);
% dataPlus = data(1,1:120000)+data(2,1:120000);
% plot(dataPlus)
% legend('EOGverPlus')
% minus = subplot(4,1,4);
% dataMinus = data(1,1:120000)-data(2,1:120000);
% plot(dataMinus)
% legend('EOGverMinus')
% linkaxes([up,down,plus,minus],'xy')
% Frontal = Fz
% Parietal = Pz
% Occipital = O1-O2
% EOGvertical = EOGup-EOGdown
chan6(1) = {'EOGvertical'};
data6(1,:) = data(strcmp('EOGup',channel),:)-data(strcmp('EOGdown',channel),:);
% EOGhorizontal = EOGright-EOGleft
chan6(2) = {'EOGhorizontal'};
data6(2,:) = data(strcmp('EOGright',channel),:)-data(strcmp('EOGleft',channel),:);
% EMG = EMG1-EMG2;
chan6(3) = {'EMG'};
data6(3,:) = data(strcmp('EMG1',channel),:)-data(strcmp('EMG2',channel),:);
% Frontal = Fz
chan6(4) = {'Frontal'};
data6(4,:) = data(strcmp('Fz',channel),:);
% Parietal = Pz
chan6(5) = {'Parietal'};
data6(5,:) = data(strcmp('Pz',channel),:);
% Occipital = O1-O2
chan6(6) = {'Occipital'};
data6(6,:) = data(strcmp('O1',channel),:)-data(strcmp('O2',channel),:);
% Frontal = data((strcmp('Fz',channel)),:);
%% Read data and segment into specified length 
% 25-second epoch - For comparison with human performance)
% 5-second epoch - For detection substages
interval = 30;
timeSeriesData = read_edf_segment(data6,fs,interval,length(chan6));% Change function name

% Number of time segment/epoch
n_seg = recordtime; % Record time (second)
n_ts = length(timeSeriesData);
%% Time segment name
for i=1:n_seg
    name = sprintf('timeseg_%d',i)
    timelabel{i}=name;
end

[labels,keywords] = labelgen(n_ts,2,chan6,timelabel);

%% Save in HCTSA input format
hctsafile = sprintf(['TS_',edffile(1:length(edffile)-9),'_',num2str(interval),'sec','_','6chan']);
% mkdir 041017_KJ_N2
% cd ./041017_KJ_N2
save(hctsafile,'timeSeriesData','labels','keywords')

