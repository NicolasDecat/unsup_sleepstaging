

%% Power spectrum: use first 37 epochs of each stage

% Load EEG Time series
load('HCTSA_N.mat','TimeSeries')
load('ccshs_1800005_annot.mat','sleepstage')

% EEG channel only 
EEG_ch = 1:length(sleepstage);  

% EEG Time Series
EEG_TS = TimeSeries(EEG_ch,4);  
EEG_TS = table2array(EEG_TS);

% Average EEG signal across all epochs for each stage
Wake_TS = EEG_TS(sleepstage == 0,:);
N1_TS = EEG_TS(sleepstage == 1,:);
N2_TS = EEG_TS(sleepstage == 2,:);
N3_TS = EEG_TS(sleepstage == 3,:);
rem_TS = EEG_TS(sleepstage == 5,:);

% Add rows next to each other to reshape in 1 row
Wake_TS2 = reshape(Wake_TS.',1,[]);
N1_TS2 = reshape(N1_TS.',1,[]);
N2_TS2 = reshape(N2_TS.',1,[]);
N3_TS2 = reshape(N3_TS.',1,[]);
rem_TS2 = reshape(rem_TS.',1,[]);

% Trim each stage to length of shortest stage (N1)
% Will need to change how it was trimmed (we're chosing the first epochs
% instead of taking randomly)
MinLength = length(N1_TS2);

Wake_TS2 = Wake_TS2(1:MinLength);
N1_TS2 = N1_TS2(1:MinLength);
N2_TS2 = N2_TS2(1:MinLength);
N3_TS2 = N3_TS2(1:MinLength);
rem_TS2 = rem_TS2(1:MinLength);

SIGNAL = [{Wake_TS2} {N1_TS2} {N2_TS2} {N3_TS2} {rem_TS2}];


%%%%%%%% Power spectrum

DecibelsFlag = 1;   % power will be converted to db

SamplingRate = 100;   

for st = 1:5
    
    % signal
    signal = SIGNAL{st};  

    % get power
    pow = (abs(fft(signal)).^2)/length(signal);

    % convert to decibels
    if DecibelsFlag==1
        pow = 10*log10(pow/max(pow));
    end
    
%     % convert to relative power (%)
%     pwrTot = bandpower(pow,SamplingRate,[0 SamplingRate/2]);
%     pow = pow/pwrTot*100;

    % first half of data without negative frequencies
    power(st,:) = pow(1:min(floor(length(signal)/2)+1,length(signal)));

    % define df and fNQ
    df = 1/(length(signal)/SamplingRate);
    fNQ = SamplingRate/2;

    faxis(st,:) = (0:df:fNQ);

end

%%% Smoothing: average every 10 frequency to make it more readable

%%% Power
for st = 1:5
    x = power(st,:)';
    S = numel(x);
    xx = reshape(x(1:S - mod(S, 100)), 100, []);
    y(:,st)  = sum(xx, 1).' / 100;
end
y = y';  
power = y;

%%% Freq
for st = 1:5
    a = faxis(st,:)';
    B = numel(a);
    aa = reshape(a(1:B - mod(S, 100)), 100, []);
    k(:,st)  = sum(aa, 1).' / 100;
end
k = k';  
faxis = k;

row10 = 1:10:length(x);
row10(end) = [];


%%% Plot power spectrum

figure;
for st = 1:5
    plot(faxis(st,:), power(st,:));
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
    hold on
end

legend('Wake','N1','N2','N3','REM')
% ylim([-120 -30])  % Random


%% Using eeglab


% Load EEG Time series
load('HCTSA_N.mat','TimeSeries')
load('ccshs_1800005_annot.mat','sleepstage')

% EEG channel only 
EEG_ch = 1:length(sleepstage);  

% EEG Time Series
EEG_TS = TimeSeries(EEG_ch,4);  
EEG_TS = table2array(EEG_TS);

% Average EEG signal across all epochs for each stage
Wake_TS = EEG_TS(sleepstage == 0,:);
N1_TS = EEG_TS(sleepstage == 1,:);
N2_TS = EEG_TS(sleepstage == 2,:);
N3_TS = EEG_TS(sleepstage == 3,:);
rem_TS = EEG_TS(sleepstage == 5,:);

% Add rows next to each other to reshape in 1 row
Wake_TS2 = reshape(Wake_TS.',1,[]);
N1_TS2 = reshape(N1_TS.',1,[]);
N2_TS2 = reshape(N2_TS.',1,[]);
N3_TS2 = reshape(N3_TS.',1,[]);
rem_TS2 = reshape(rem_TS.',1,[]);

% Trim each stage to length of shortest stage (N1)
% Will need to change how it was trimmed (we're chosing the first epochs
% instead of taking randomly)
MinLength = length(N1_TS2);

Wake_TS2 = Wake_TS2(1:MinLength);
N1_TS2 = N1_TS2(1:MinLength);
N2_TS2 = N2_TS2(1:MinLength);
N3_TS2 = N3_TS2(1:MinLength);
rem_TS2 = rem_TS2(1:MinLength);

% % Calculate spectral power - 1 stage
% figure;
% [spectra,freqs] = spectopo(Wake_TS2, 0, 1000);

% Calculate spectral power - all stages
figure;
ALL = [Wake_TS2 ; N1_TS2 ; N2_TS2 ; N3_TS2 ; rem_TS2];
[spectra,freqs] = spectopo(ALL, 0, 1000);
legend('Wake','N1','N2','N3','rem');


%% Get misclassified N2 epochs and plot their spectral power

% First, run hctsa_allfeatures (1 sub, 1 chan, 1 Nf)
%run('/Users/nico/Documents/GitHub/unsup_sleepstaging/hctsa_scripts/hctsa_allfeatures')

% Load TimeSeries
load('HCTSA_N.mat','TimeSeries')
Data = TimeSeries{:,4};

%%%% Original labels N2
orig_labels_idx = statsOut.scoredTest(1,23:33);

%%%% Labels given by algo to N2 epochs
algo_labels_idx = statsOut.predictTest(1,23:33);

% Testing epochs labelled N2 by original labels
N2_epoch_ID = testTS(1,23:33); 
N2_epoch_EEG = Data(N2_epoch_ID,:);  % 11 testing epochs 

% Get EEG of epochs not labelled correctly by algo (not labelled N2)
idx_not_N2 = find(algo_labels_idx ~= 2);
EEG_Not_N2 = N2_epoch_EEG(idx_not_N2,:);

% Get EEG of epochs  labelled correctly by algo ( labelled N2)
idx_N2 = find(algo_labels_idx == 2);
EEG_N2 = N2_epoch_EEG(idx_N2,:);

% Plot spectral power of (mis)classified N2 epochs
PlotWrongN2 = true;

if PlotWrongN2
    for i = 1:size(EEG_Not_N2,1)
        figure;
        [spectra,freqs] = spectopo(EEG_Not_N2(i,:), 0, 1000);
        legend('Misclassified N2 epoch');
    end
else
    for i = 1:size(EEG_N2,1)
        figure;
        [spectra,freqs] = spectopo(EEG_N2(i,:), 0, 1000);
        legend('Misclassified N2 epoch');
    end
end




