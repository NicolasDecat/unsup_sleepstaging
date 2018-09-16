clc;

% %% CONFIGURATION
EDF_FILE='/Volumes/Spaceship/TS_ccshs-trec-1800749_7chan.mat';
ANSWER_FILE = '/Volumes/Spaceship/ccshs_datasets/ccshs_1800749_annot.mat';
% Global

load(EDF_FILE);

%% Extract channel names, sampling frequency from header struct
fs = 128;

ks = split(string(keywords)', ',');
ksu = unique(ks(:,1), 'stable');
nchannels = length(ksu);
data = timeSeriesData;

one_channel_length = size(data, 1)/nchannels;
old_data = data;
data=reshape(old_data', one_channel_length*size(old_data,2), nchannels)';

channel_info=[(1:length(ksu))',string(ksu)]

channels.EEG = [1];
channels.EOG = [2];
channels.EMG = [3];
    
chans=channels;
data = data.* 1000;
%%
% samplingrate=200; %sampling rate in Hertz
% mspersample=1000/samplingrate; %milliseconds per sample 
% eodata=data([1 9],:);
% eodata=eodata';
% mytime=[1:1:size(eodata,1)]'; %create column vector of length same as your data 
% mytime=mytime*mspersample; %adjust time vector to sample rate
% data= [mytime eodata]; %attach the time column to your x and y data already stored in 'data'

annotation = load(ANSWER_FILE);
label = annotation.sleepstage;

assert(size(label, 1) == one_channel_length, "The length of the annotation must be the same as one channel length");
%% Pre-processing

%% Downsample
% ds=[];
% for i=1:size(data,1)
%    ds=[ds; downsample(data(i, :), 2)];
% end
% data=ds;
% data = medfilt1(data,3);
% fs = fs/2;
es = 30;

%% Analyse each epochs


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
    draw_subplot(ax, data, "EEG", EEG_CHN, cellstr(channel_info(EEG_CHN, 2))', no_of_epochs, [-300, 300], fs, es, chans);

    ax=subplot(3, 1, 2);
    draw_subplot(ax, data, "EOG", EOG_CHN, cellstr(channel_info(EOG_CHN, 2))', no_of_epochs, [-300, 300], fs, es, chans);

    ax=subplot(3, 1, 3);
    draw_subplot(ax, data, "EMG", EMG_CHN, cellstr(channel_info(EMG_CHN, 2))', no_of_epochs, [-300, 300], fs, es, chans);

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
        
        plt=plot(ax, frame*(no_of_epochs-1)+1:frame*no_of_epochs, ts_data(:, frame*(no_of_epochs-1)+1:frame*no_of_epochs));

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
            ax.YAxisLocation = 'origin';
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