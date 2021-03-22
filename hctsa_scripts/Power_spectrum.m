
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

b = figure;
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
WakeC = {cb(1,:)};
N1C = {cb(12,:)};
N2C = {cb(3,:)};
N3C = {cb(4,:)};
REMC = {cb(5,:)};
Color = [WakeC N1C N2C N3C REMC];

% Reorder: N1 N2 N3 Wake REM
SIGNALS(1) = SIGNAL(2);
SIGNALS(2) = SIGNAL(3);
SIGNALS(3) = SIGNAL(1);
SIGNALS(4) = SIGNAL(4);
SIGNALS(5) = SIGNAL(5);

for i = 1:5
    signal = SIGNALS{i};
    pwelch(signal,w_window,w_overlap,freqV,Fs,'psd');
    h = get(gca, 'Children');
    set(h(1), 'Color', Color{i});
    set(h(1), 'LineWidth', 1.3);
    hold on
     
end

ylabel('Power (dB)')
title('Power spectrum of N2 epochs (as labelled by experts)')
legend('Wake','N1','N2','N3','REM')

% set(b, 'Color', 'w')
% fpath = '/Users/nico/Documents/HCTSA/Analysis/spectral';
% export_fig([fpath filesep 'spectrum_N2_epochs'],'-r 300')

%%%% Bar graph for each frequency band
load('/Users/nico/Documents/HCTSA/Analysis/spectral/signal.mat')
Fs = 128;
w_window=40*Fs; % Where Fs is sampling rate (applies the method on windows of 6s)
w_overlap=w_window/2; % overlap between windows
df=0.2; % freq resolution
freqV=1:0.2:40;

freq_band = {'delta (1-4 Hz)','theta (4-8 Hz)','sigma (11-16 Hz)','beta (20-40 Hz)'};
FREQ = {[1 4] [4 8] [11 16] [20 40]};

load('/Users/nico/Documents/HCTSA/Analysis/spectral/signal.mat')

figure;
[ha, pos] = tight_subplot(2,4,[.1 .1],[.1 .05],[.05 .05]);


for i = 1:4   % for each frequency
    
    axes(ha(i+4)); 
    
    for j = 1:5   % for each stage
       
        % Get the epochs from one stage and compute their power
        signal = SIGNALS{j};
        [pow,f] = pwelch(signal,w_window,w_overlap,freqV,Fs,'psd');

        % Idx of data points within the frequency band
        within_band = find(f >= FREQ{i}(1) & f <= FREQ{i}(2));
        
        % Mean power
        Power(i,j) = mean(pow(within_band));   % Table: each row a freq, each col a stage
        
    end
      
    
    X = categorical({'Wake','N1','N2','N3','REM'});
    X = reordercats(X,{'Wake','N1','N2','N3','REM'});

    B = bar(X,Power(i,:));
    B.FaceColor = 'flat';
    B.CData = [Color{1};Color{2};Color{3};Color{4};Color{5}];
    ylabel('Power/frequency (dB/Hz)')
    title(sprintf('%s band',string(freq_band(i))))

    
     axes(ha(i));
     
     for s = 1:5
        signal = SIGNALS{s};
        pwelch(signal,w_window,w_overlap,freqV,Fs,'psd');
        h = get(gca, 'Children');
        set(h(1), 'Color', Color{s});
        set(h(1), 'LineWidth', 1.3);
        hold on
     end
     
     title(sprintf('Welch PSD (%s)',string(freq_band(i))))
     ylabel('Power (dB)')
     xlim([FREQ{i}(1) FREQ{i}(2)])
     
    % Legend to first plot
    if i == 1
        legend('Wake','N1','N2','N3','REM')
    end
    
end

%%% Same figures but normalize power across bands

for i = 1:5   % for each stage, get total power across the 4 bands
    ttl_pow(1,i) =  sum(Power(:,i)); 
    power_norm(:,i) = Power(:,i)./ttl_pow(i);
end

a = figure;
[ha, pos] = tight_subplot(2,4,[.1 .1],[.1 .05],[.05 .05]);


for i = 1:4   % for each frequency
    
    axes(ha(i+4)); 
    
    for j = 1:5   % for each stage
       
        % Get the epochs from one stage and compute their power
        signal = SIGNALS{j};
        [pow,f] = pwelch(signal,w_window,w_overlap,freqV,Fs,'psd');

        % Idx of data points within the frequency band
        within_band = find(f >= FREQ{i}(1) & f <= FREQ{i}(2));
        
        % Mean power
        Power(i,j) = mean(pow(within_band));   % Table: each row a freq, each col a stage
        
    end
      
    
    X = categorical({'Wake','N1','N2','N3','REM'});
    X = reordercats(X,{'Wake','N1','N2','N3','REM'});

    B = bar(X,power_norm(i,:));    
    B.FaceColor = 'flat';
    B.CData = [Color{1};Color{2};Color{3};Color{4};Color{5}];
    ylabel('Relative power (%)')
    title(sprintf('%s band',string(freq_band(i))))

    
     axes(ha(i));
     
     for s = 1:5
        signal = SIGNALS{s};
        pwelch(signal,w_window,w_overlap,freqV,Fs,'psd');
        h = get(gca, 'Children');
        set(h(1), 'Color', Color{s});
        set(h(1), 'LineWidth', 1.3);
        hold on
     end
     
     title(sprintf('Welch PSD (%s)',string(freq_band(i))))
     ylabel('Power (dB)')
     xlim([FREQ{i}(1) FREQ{i}(2)])
     
    % Legend to first plot
    if i == 1
        legend('Wake','N1','N2','N3','REM')
    end
    
end

% set(a, 'Color', 'w')
% fpath = '/Users/nico/Documents/HCTSA/Analysis/spectral';
% export_fig([fpath filesep 'power_per_band_relative'],'-r 300')




%% Using Welch (in db instead of db/hz)

load('/Users/nico/Documents/HCTSA/Analysis/spectral/signal.mat')


Fs = 128;
w_window=40*Fs; % Where Fs is sampling rate (applies the method on windows of 6s)
w_overlap=w_window/2; % overlap between windows
df=0.2; % freq resolution
freqV=1:0.2:40;

figure;
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
WakeC = {cb(1,:)};
N1C = {cb(12,:)};
N2C = {cb(3,:)};
N3C = {cb(4,:)};
REMC = {cb(5,:)};
Color = [N2C WakeC N1C  N3C REMC];

figure; 
for i = 1:5
    signal = SIGNAL{i};
    [pow(i,:),f(i,:)] = pwelch(signal,w_window,w_overlap,freqV,Fs,'psd');   
end

pow = pow2db(pow);

figure;
for i = 1:5
    plot(f(i,:),pow(i,:))
    h = get(gca, 'Children');
    set(h(1), 'Color', Color{i});
    set(h(1), 'LineWidth', 1.3);
    hold on
end

title('Power spectrum of N2 epochs (as labelled by experts)')
legend('N2','Wake','N1','N3','REM')

% Bar

freq_band = {'delta (1-4 Hz)','theta (4-8 Hz)','alpha (8-16 Hz)','beta (20-40 Hz)'};
FREQ = {[1 4] [4 8] [8 16] [20 40]};

figure; ax = gca;
[ha, pos] = tight_subplot(1,4,[.045 .1],[.1 .05],[.05 .05]);


for i = 1:4   % for each frequency
    
    axes(ha(i)); 
    
    for j = 1:5   % for each stage
       
        % Idx of data points within the frequency band
        within_band = find(f(j,:) >= FREQ{i}(1) & f(j,:) <= FREQ{i}(2));
        
        % Mean power
        Power(i,j) = mean(pow(j,within_band));   % Table: each row a freq, each col a stage
     
    end
      
    X = categorical({'N2','Wake','N1','N3','REM'});
    X = reordercats(X,{'N2','Wake','N1','N3','REM'});

    B = bar(X,Power(i,:));
    B.FaceColor = 'flat';
    B.CData = [Color{1};Color{2};Color{3};Color{4};Color{5}];
    ylabel('Power (dB)')
    set(gca, 'YDir','reverse')
    title(sprintf('%s Power',string(freq_band(i))))

end
    







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
    Power(st,:) = pow(1:min(floor(length(signal)/2)+1,length(signal)));

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
        x = Power(st,:)';
        S = numel(x);
        xx = reshape(x(1:S - mod(S, 128)), 128, []);
        y(:,st)  = sum(xx, 1).' / 128;
    end
    y = y';  
    Power = y;

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
    
    plot(faxis(st,:), Power(st,:));
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

