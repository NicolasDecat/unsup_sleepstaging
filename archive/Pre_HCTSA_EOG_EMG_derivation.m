%% Preprocessing of polysomnographic data and generate time-series matrix for HCTSA
% Before submitting data for features extraction,
% 1. Read data in EDF format into matrix/cell array of 30-second epochs
% 2. Generate labels and keywords to identify time-series
% 3. Save data in the format suitable for HCTSA
% #########################################################################
%
configuration_settings;

whichData = WHICH_DATA;
edfname = EDF_FILE;
edffile = strcat(DATA_DIR, filesep, EDF_FILE);
%% Specify case (How many channel?)
nChannel = NUM_CHANNELS; % Number of channels to be used:[1,2,3] 
% Read data and segment into specified length 
% 25-second epoch - For comparison with human performance)
% 5-second epoch - For detection substages
interval = 30; % unit: second
%% Load data using blockEdfLoad
% https://sleepdata.org/tools/dennisdean-block-edf-loader 
addpath(genpath(BLOCKEDFLOAD_DIR));
[header, signalHeader, signalCell] = blockEdfLoad(edffile);

%% Extract channel names, sampling frequency for different cases

C3 = 1;
C4 = 2;
A1 = 3;
A2 = 4;
LOC = 5;
ROC = 6;
EMG1 = 13;
EMG2 = 14;

% EOG
% selectedSignal(1).raw = signalCell{LOC}-signalCell{A2};
% selectedHeader(1).signal_labels = strcat(signalHeader(LOC).signal_labels,'-',signalHeader(A2).signal_labels);
% selectedHeader(1).sampling_rate = signalHeader(LOC).samples_in_record;
% 
% selectedSignal(2).raw = signalCell{ROC}-signalCell{A1};
% selectedHeader(2).signal_labels = strcat(signalHeader(ROC).signal_labels,'-',signalHeader(A1).signal_labels);
% selectedHeader(2).sampling_rate = signalHeader(ROC).samples_in_record;
% 
% selectedSignal(3).raw = (signalCell{LOC}-signalCell{A2}) + (signalCell{ROC}-signalCell{A1});
% selectedHeader(3).signal_labels = strcat(signalHeader(LOC).signal_labels,'-',signalHeader(A2).signal_labels, '+', signalHeader(ROC).signal_labels,'-',signalHeader(A1).signal_labels);
% selectedHeader(3).sampling_rate = signalHeader(LOC).samples_in_record;
% 
% selectedSignal(4).raw = (signalCell{LOC}-signalCell{A2}) - (signalCell{ROC}-signalCell{A1});
% selectedHeader(4).signal_labels = strcat(signalHeader(LOC).signal_labels,'-',signalHeader(A2).signal_labels, '-', signalHeader(ROC).signal_labels,'-',signalHeader(A1).signal_labels);
% selectedHeader(4).sampling_rate = signalHeader(LOC).samples_in_record;
% 
% selectedSignal(5).raw = (signalCell{LOC}-signalCell{A2}) .* (signalCell{ROC}-signalCell{A1});
% selectedHeader(5).signal_labels = strcat(signalHeader(LOC).signal_labels,'-',signalHeader(A2).signal_labels, '*', signalHeader(ROC).signal_labels,'-',signalHeader(A1).signal_labels);
% selectedHeader(5).sampling_rate = signalHeader(LOC).samples_in_record;

% EMG
selectedSignal(1).raw = signalCell{EMG1};
selectedHeader(1).signal_labels = strcat(signalHeader(EMG1).signal_labels);
selectedHeader(1).sampling_rate = signalHeader(EMG1).samples_in_record;

selectedSignal(2).raw = signalCell{EMG2};
selectedHeader(2).signal_labels = strcat(signalHeader(EMG2).signal_labels);
selectedHeader(2).sampling_rate = signalHeader(EMG2).samples_in_record;

selectedSignal(3).raw = signalCell{EMG1} + signalCell{EMG2};
selectedHeader(3).signal_labels = strcat(signalHeader(EMG1).signal_labels, '+', signalHeader(EMG2).signal_labels);
selectedHeader(3).sampling_rate = signalHeader(EMG1).samples_in_record;

selectedSignal(4).raw = signalCell{EMG1} - signalCell{EMG2};
selectedHeader(4).signal_labels = strcat(signalHeader(EMG1).signal_labels, '-', signalHeader(EMG2).signal_labels);
selectedHeader(4).sampling_rate = signalHeader(EMG1).samples_in_record;

selectedSignal(5).raw = signalCell{EMG1} .* signalCell{EMG2};
selectedHeader(5).signal_labels = strcat(signalHeader(EMG1).signal_labels, '.*', signalHeader(EMG2).signal_labels);
selectedHeader(5).sampling_rate = signalHeader(EMG1).samples_in_record;


% Downsample EMG channel
% selectedSignal(3).raw = downsample(selectedSignal(3).raw,2); % Downsampled by factor of 2
% selectedHeader(3).sampling_rate = selectedHeader(3).sampling_rate/2;

nChannel = length(selectedSignal);

% Segmentation
for m = 1:nChannel
    selectedSignal(m).chopped = read_edf_segment(selectedSignal(m).raw',...
        selectedHeader(m).sampling_rate,interval,1);
    % Generate time labels
    [n_ts,~]  = size(selectedSignal(m).chopped);
    for t = 1:n_ts
        name = sprintf('timeseg_%d',t);
        timelabel{t} = name;
    end

    % Combine 2/3 channels into single timeSeriesData matrix
    timeSeriesData(n_ts*(m-1)+1:n_ts*m,:) = selectedSignal(m).chopped;
    [labels(n_ts*(m-1)+1:n_ts*m),keywords(n_ts*(m-1)+1:n_ts*m)] = labelgen(n_ts,2,{selectedHeader(m).signal_labels},timelabel);            
end
        
%% Save in HCTSA input format
hctsafile = strcat('TS_',edfname(1:length(edfname)-4),'_EMG_',num2str(nChannel),'chan');
save(hctsafile,'timeSeriesData','labels','keywords')

