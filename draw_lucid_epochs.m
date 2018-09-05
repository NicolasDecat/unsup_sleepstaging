clc;

% %% CONFIGURATION
% EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/KJ_N1_data.mat';
% EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/ME_N1_data.mat';
% EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/KJ_N2_data.mat';
%EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/ME_N1_data.mat';
EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/ME_N2_data.mat';
% EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/ME_N3_data.mat';
% EDF_FILE='/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/LI_N2_data.mat';
%EDF_FILE='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging_master/Voss_Pilot/SK_N1_F1_data.mat';
%EDF_FILE='/Users/Zhao/SleepPsychoPhysics/Source/unsup_sleepstaging_master/Voss_Pilot/SK_N1_F1_data.mat';


% Global

load(EDF_FILE);

%% Extract channel names, sampling frequency from header struct
channel = cellstr(header.channelname);  % Channel names + converted into cell array
fs = header.samplerate(1);      % Assume sampling rate are the same for all channels
recordtime = header.records;    % aaTotal recorded time (seconds)
nchannels = header.channels;    % Number of channels

% Amplitude calibration
phys_range = header.physmax - header.physmin;
dig_range = header.digimax - header.digimin;

gain = phys_range ./ dig_range;

for i = 1:size(data, 1)
   data(i, :) = (data(i, :) - header.digimin(i)) * gain(i) + header.physmin(i);
end

channel_info=[(1:length(header.channelname))',string(header.channelname)]

%% 
% KJ_N1
dataset_names=["KJ_N1", "KJ_N2", "LI_N2", "ME_N1", "ME_N2", "ME_N3"];
% dataset=find(strcmp(dataset_names,"KJ_N1"));
% data(6,:)=data(5,:)-data(6,:);
% data(3,:)=data(3,:)-data(24,:);
% data(4,:)=data(4,:)-data(24,:);

channels{find(strcmp(dataset_names,"KJ_N1"))}.EEG = [11, 12];
channels{find(strcmp(dataset_names,"KJ_N1"))}.EOG = [9,10];
channels{find(strcmp(dataset_names,"KJ_N1"))}.EMG = [6];
    
% KJ_N2
% dataset=find(strcmp(dataset_names,"KJ_N2"));
% data(9,:)=data(9,:)-data(24,:);
% data(10,:)=data(10,:)-data(25,:);
% data(3,:)=data(3,:)-data(4,:);
channels{find(strcmp(dataset_names,"KJ_N2"))}.EEG = [9,10];
channels{find(strcmp(dataset_names,"KJ_N2"))}.EOG = [7,8];
channels{find(strcmp(dataset_names,"KJ_N2"))}.EMG = [3];

% % LI_N2
% dataset=find(strcmp(dataset_names,"LI_N2"));
% data(3,:)=data(3,:)-data(4,:);
channels{find(strcmp(dataset_names,"LI_N2"))}.EEG = [9,10];
channels{find(strcmp(dataset_names,"LI_N2"))}.EOG = [1,2];
channels{find(strcmp(dataset_names,"LI_N2"))}.EMG = [3];


% % ME_N1
% dataset=find(strcmp(dataset_names,"ME_N1"));
% data(3,:)=data(3,:)-data(4,:);
% data(9,:)=data(9,:)-data(24,:);
% data(10,:)=data(10,:)-data(25,:);
channels{find(strcmp(dataset_names,"ME_N1"))}.EEG = [9,10];
channels{find(strcmp(dataset_names,"ME_N1"))}.EOG = [1,2];
channels{find(strcmp(dataset_names,"ME_N1"))}.EMG = [3];

% % ME_N2
dataset=find(strcmp(dataset_names,"ME_N2"));
data(3,:)=data(3,:)-data(4,:);
channels{find(strcmp(dataset_names,"ME_N2"))}.EEG = [7,8];
channels{find(strcmp(dataset_names,"ME_N2"))}.EOG = [1,2];
channels{find(strcmp(dataset_names,"ME_N2"))}.EMG = [3];

% % ME_N3
% dataset=find(strcmp(dataset_names,"ME_N3"));
% data(3,:)=data(3,:)-data(4,:);
channels{find(strcmp(dataset_names,"ME_N3"))}.EEG = [9,10];
channels{find(strcmp(dataset_names,"ME_N3"))}.EOG = [1,2];
channels{find(strcmp(dataset_names,"ME_N3"))}.EMG = [3];


chans=channels{dataset};

%% Draw power spectrum
% openExample('signal/FrequencyAnalysisExample');
% data_test=data(11,:)';
% pxx = pwelch(data_test);
% pxx = 10*log10(pxx);
% 
% NFFT=length(data_test);
% Y = fft(data_test, NFFT);
% F = ((0:1/NFFT:1-1/NFFT)*fs).';
% 
% magnitudeY = abs(Y);        % Magnitude of the FFT
% phaseY = unwrap(angle(Y));  % Phase of the FFT
% 
% helperFrequencyAnalysisPlot1(F,magnitudeY,phaseY,NFFT);

%%
% samplingrate=200; %sampling rate in Hertz
% mspersample=1000/samplingrate; %milliseconds per sample 
% eodata=data([1 9],:);
% eodata=eodata';
% mytime=[1:1:size(eodata,1)]'; %create column vector of length same as your data 
% mytime=mytime*mspersample; %adjust time vector to sample rate
% data= [mytime eodata]; %attach the time column to your x and y data already stored in 'data'

%% Pre-processing

% Notch-filter of 50Hz

% Design a filter with a Q-factor of Q=35 to remove a 50 Hz tone from 
% system running at 300 Hz.
Wo = 50/(fs/2);  BW = Wo/35;
[b,a] = iirnotch(Wo,BW);
data(chans.EEG,:)=filter(b,a,data(chans.EEG,:));

% Bandpass between [0.5 70]
[b,a]=butter(2,[0.5,70]/(fs/2),'bandpass');
sub_ts=data(chans.EEG,:)';
sub_ts=filtfilt(b,a,sub_ts);
data(chans.EEG,:) = sub_ts';

% [b,a]=butter(2,[0.5,15]/(fs/2),'bandpass');
% sub_ts=data(chans.EOG,:)';
% sub_ts=filtfilt(b,a,sub_ts);
% data(chans.EOG,:) = sub_ts';

% [b,a]=butter(6,[0.5,30]/(fs/2),'bandpass');
% sub_ts=data(chans.EOG,:)';
% sub_ts=filtfilt(b,a,sub_ts);
% data(chans.EOG,:) = sub_ts';


%%
ds=[];
for i=1:size(data,1)
   ds=[ds; downsample(data(i, :), 2)];
end
data=ds;
% data = medfilt1(data,3);
fs = fs/2;
es = 30;

%% Analyse each epochs

% Analyse consecutive high amplitude

last_for_seconds=2;
eog_max_threshold=150;
eog_boundlimit = 200;
eog_min_threshold=-150;

low_eog_max_threshold=50;
low_eog_min_threshold=-50;
consecutive_low = fs*0.8;

low_emg_max_threshold=5.5;
low_emg_min_threshold=-5.5;
consecutive_emg_per_unit=fs*0.95;
total_percentage_low_emg_per_epoch=0.90;

seconds_splitting=1;

%% Prepare for the data to be reshape

trailing_epochs = size(data, 2) - (floor(size(data, 2)/(fs*es)))*(fs*es);
data=data(:, 1:(size(data, 2)-trailing_epochs));

%%
eog_epochs=[];
% for i=1:length(chans.EOG)
    eog_chan_id=chans.EOG(1);
%     fprintf('\nAnalysing channel %s...\n', channel_info(chan_id, 2));
    
    timeframe = reshape(data(eog_chan_id, :)', es*fs, floor(size(data(eog_chan_id, :), 2)/es/fs))';
    
    [epochs, total_unit_value] = size(timeframe);
    for j=1:epochs
        epoch = timeframe(j, :);
        time_matrix = reshape(epoch, fs/seconds_splitting, es*seconds_splitting)';
%         time_matrix(find(time_matrix > eog_boundlimit)) = eog_boundlimit; 
%         time_matrix(find(time_matrix < eog_boundlimit * -1)) = eog_boundlimit*-1; 
        
        %mean_time = mean(abs(time_matrix)')
%         idx = [true,mean_time(1:end-1)>eog_max_threshold] & mean_time>eog_max_threshold & [mean_time(2:end)>eog_max_threshold,true];
%         idx = [false,idx(1:end-1)] | idx | [idx(2:end),false];
%         
%         idx2 = [true,mean_time(1:end-1)>eog_min_threshold] & mean_time>eog_min_threshold & [mean_time(2:end)>eog_min_threshold,true];
%         idx2 = [false,idx2(1:end-1)] | idx2 | [idx2(2:end),false];

        tmat = time_matrix > eog_max_threshold | time_matrix < eog_min_threshold;
        idx = [ones(size(tmat, 1), 1), tmat(:,1:end-1) == 1] & (tmat == 1) & [tmat(:,2:end) == 1, ones(size(tmat, 1), 1)];
        
        for k=1:size(idx,1)
           idx(k,:) = [false,idx(k,1:end-1)] | idx(k,:) | [idx(k,2:end),false];
        end
        
        s=sum(idx');
        idx2 = [true, s(1:end-1) > 0] & (s>0) & [s(2:end) > 0, true];
        idx2 = [false, idx2(1:end-1)] | idx2 | [idx2(2:end), false];
%         
        % Calculate consecutive low EOG signal
        tmat_low = time_matrix < low_eog_max_threshold & time_matrix > low_eog_min_threshold;
        idx = [ones(size(tmat_low, 1), 1), tmat_low(:,1:end-1) == 1] & (tmat_low == 1) & [tmat_low(:,2:end) == 1, ones(size(tmat_low, 1), 1)];
        
        for k=1:size(idx,1)
           idx(k,:) = [false,idx(k,1:end-1)] | idx(k,:) | [idx(k,2:end),false];
        end
        
        ss=sum(idx');
        idx3 = [true, ss(1:end-1) > consecutive_low] & (ss> consecutive_low) & [ss(2:end) > consecutive_low, true];
        idx3 = [false, idx3(1:end-1)] | idx3 | [idx3(2:end), false];
        
        if sum(idx2) > last_for_seconds && (length(find(idx3==1))/length(idx3)) > 0.5
            eog_epochs = [eog_epochs j];
        end
    end
    
    emg_chan_id=chans.EMG(1);
    timeframe = reshape(data(emg_chan_id, :)', es*fs, size(data(emg_chan_id, :), 2)/es/fs)';

    final_epochs = [];
    for i = 1:length(eog_epochs)
        epoch = timeframe(eog_epochs(i), :);
        time_matrix = reshape(epoch, fs/seconds_splitting, es*seconds_splitting)';
        tmat=time_matrix < low_emg_max_threshold & time_matrix > low_emg_min_threshold;
        tmats=sum(tmat');
        
        consecutive_low_emg = length(find(tmats >= consecutive_emg_per_unit));
        if (consecutive_low_emg/size(tmats,2) >= total_percentage_low_emg_per_epoch)
            final_epochs = [final_epochs eog_epochs(i)];
        end
    end
    
    for k=1:length(final_epochs)
        fprintf('Epoch %d starting with %d seconds.\n', final_epochs(k), ...
            ((final_epochs(k)-1) * es));
    end
%end



%% Plot the relevant channel
f=figure('rend','painters','pos',[10 10 1200 800]);
set(f,'color','w');
draw(f, channel_info, data, 0, fs, es, chans);

% Create slider
sld = uicontrol('Parent',f,'Style', 'slider',...
    'Min',0,'Max',size(data,2)/fs,'Value',0,...
    'Position', [320 0.05 600 35],...
    'FontSize', 10, ...
    'Visible', 'on', ...
    'SliderStep', [1/((size(data,2)/fs)/es), 0.1],...
    'Callback', {@sliderselected, channel_info, data, fs, es, chans});
tb=uicontrol('Parent',f,'Style', 'edit', 'Position', [240 25 50 25],'Callback', {@tbentered, channel_info, data, sld, fs, es, chans});

warning('off');

function draw(f, channel_info, data, num_seconds, fs, es, chans)
    no_of_epochs=(num_seconds/es)+1;
    EEG_CHN = chans.EEG;
    EOG_CHN = chans.EOG;
    EMG_CHN = chans.EMG;

    ax=subplot(3, 1, 1);
    draw_subplot(ax, data, "EEG", EEG_CHN, cellstr(channel_info(EEG_CHN, 2))', no_of_epochs, [-200 200], fs, es, chans);

    ax=subplot(3, 1, 2);
    draw_subplot(ax, data, "EOG", EOG_CHN, cellstr(channel_info(EOG_CHN, 2))', no_of_epochs, [-200, 200], fs, es, chans);

    ax=subplot(3, 1, 3);
    draw_subplot(ax, data, "EMG", EMG_CHN, cellstr(channel_info(EMG_CHN, 2))', no_of_epochs, [-10, 10], fs, es, chans);

    ButtonP=uicontrol('Parent',f,'Style','pushbutton','String','Previous','Units','normalized','Position',[0.05 0.02 0.1 0.05],'Visible','on','Callback',{@plotbuttonclick, channel_info, data, fs, es, chans, num_seconds-es});
    ButtonN=uicontrol('Parent',f,'Style','pushbutton','String','Next','Units','normalized','Position',[0.80 0.02 0.1 0.05],'Visible','on','Callback',{@plotbuttonclick, channel_info, data, fs, es, chans, num_seconds+es});
end

function draw_subplot(ax, data, plot_title, chnls, channel_labels, no_of_epochs, y_limits, fs, es, chans)
    pos=get(ax, 'position');
    pos(4) = 0.26;
    set(ax, 'position', pos);
    
    for i=1:length(chnls)
        chn_idx=chnls(i);
        ts_data = data(chn_idx,:);

        frame=es*fs;
        
        if chn_idx==15 || chn_idx==16
            plt=plot(ax, frame*(no_of_epochs-1)+1:frame*no_of_epochs, ts_data(:, frame*(no_of_epochs-1)+1:frame*no_of_epochs), 'LineWidth', 0.5);
        else
            plt=plot(ax, frame*(no_of_epochs-1)+1:frame*no_of_epochs, ts_data(:, frame*(no_of_epochs-1)+1:frame*no_of_epochs), 'LineWidth', 0.5);
        end
        hold on;

        axis tight;
        ax.YGrid = 'on';
        ax.XGrid = 'on';

        pos = get(ax, 'Position');
        pos(1) = 0.055;
        pos(3) = 0.9;
        set(ax, 'Position', pos)

        set(ax,'box','off');
        ylabel(ax, plot_title);
        legend(channel_labels);
        legend(ax, 'show');
        set(ax, 'XTick', 0:fs:size(data, 2));
        XTickLabel = get(ax,'XTick');
        set(ax, 'XTickLabel', string(XTickLabel/fs));
        
%         YTickLabel = get(ax,'YTick');
%         set(ax, 'YTickLabel', string(YTickLabel));

        if (length(y_limits) > 0)
            %ax.YAxisLocation = 'origin';
            ylim(ax, y_limits);
        end
    end
end

function plotbuttonclick(h,~,channel_info, data, fs, es, chans, num_seconds)
    f = h.Parent;
    clearsubplots(f);
    draw(f, channel_info, data, num_seconds, fs, es, chans);
end

function sliderselected(h,~,channel_info, data, fs, es, chans)
    f = h.Parent;
    clearsubplots(f);
    draw(f, channel_info, data, h.Value, fs, es, chans);
end

function tbentered(h,~,channel_info, data, sd, fs, es, chans)
    f = h.Parent;
    clearsubplots(f);
    n=str2num(get(h,'String'));
    draw(f, channel_info, data, n, fs, es, chans);
    set(sd, 'Value', n);
    
end

function clearsubplots(f)
    h=get(f,'children');
    for i = 1:length(h)
        if (strcmp(class(h(i)),'matlab.graphics.axis.Axes'))
            cla(h(i));
        end
    end
end