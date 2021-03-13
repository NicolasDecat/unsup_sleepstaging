
%% Power spectra of (mis)classified N2 epochs

%%%% Power spectra of N2 epochs (as defined by original labels) that were
%%%% classified as N2 or as any other stage by our algorithm


%% Using EEGLAB (spectopo function)

%%% Load the time-series of the N2 epochs (80 epochs each)
%%% Row 1 = N2 (correctly classified N2 epochs), 2 = wake, 3 = N1, 4 = N3, 5 = REM
load('/Users/nico/Documents/HCTSA/Analysis/spectral/signal.mat')

% Convert cell array to numeric array, each row being a different stage
for st = 1:5
    signal(st,:) = SIGNAL{1,st};
end

% Plot power spectra
figure;

[spectra,freqs] = spectopo(signal, 0,128);
grid on
xlim([0 45])

legend('N2','Wake','N1','N3','REM')

title('Power spectrum of N2 epochs (as labelled by experts)')

%% Using Welch

load('/Users/nico/Documents/HCTSA/Analysis/spectral/signal.mat')


Fs = 128;
w_window=40*Fs; % Where Fs is sampling rate (applies the method on windows of 6s)
w_overlap=w_window/2; % overlap between windows
df=0.2; % freq resolution
freqV=1:0.2:40;

figure;

for i = 1:5
    signal = SIGNAL{i};
    pwelch(signal,w_window,w_overlap,freqV,Fs,'psd');
    plot(faxis, pow);
    hold on
end

legend('N2','Wake','N1','N3','REM')

%% Using Thomas' script

DecibelsFlag = 1;   
SamplingRate = 128;   

%%% Load the time-series of the N2 epochs (80 epochs each)
%%% Cell 1 = N2 (correctly classified N2 epochs), 2 = wake, 3 = N1, 4 = N3, 5 = REM
load('/Users/nico/Documents/HCTSA/Analysis/spectral/signal.mat')


for st = 1:5     % For all 5 stages
    
    % signal from each epoch 
    signal = SIGNAL{st};  

    % get power
    pow = (abs(fft(signal)).^2)/length(signal);

    % convert to decibels
    if DecibelsFlag==1
        pow = 10*log10(pow/max(pow));
    end

    % first half of data without negative frequencies
    power(st,:) = pow(1:min(floor(length(signal)/2)+1,length(signal)));

    % define df and fNQ
    df = 1/(length(signal)/SamplingRate);
    fNQ = SamplingRate/2;

    faxis(st,:) = (0:df:fNQ);
    
end

%%% Smoothing: average every 100 frequency data points to make the graph more readable
smooth = false;

if smooth
    
    %%% Every 100 power values are averged 
    for st = 1:5
        x = power(st,:)';
        S = numel(x);
        xx = reshape(x(1:S - mod(S, 128)), 128, []);
        y(:,st)  = sum(xx, 1).' / 128;
    end
    y = y';  
    power = y;

    %%% Every 100 frequency values are averaged
    for st = 1:5
        a = faxis(st,:)';
        B = numel(a);
        aa = reshape(a(1:B - mod(S, 128)), 128, []);
        k(:,st)  = sum(aa, 1).' / 128;
    end
    k = k';  
    faxis = k;

    row10 = 1:10:length(x);
    row10(end) = [];
end


%%% Plot power spectrum
figure;

for st = 1:5
    
    plot(faxis(st,:), power(st,:));
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
    hold on
    
end

legend('Epochs labelled as N2 by the algorithm',...
    'Epochs labelled as Wake by the algorithm',...
    'Epochs labelled as N1 by the algorithm',...
    'Epochs labelled as N3 by the algorithm',...
    'Epochs labelled as REM by the algorithm');

title('Power spectrum of N2 epochs (as labelled by experts)')

