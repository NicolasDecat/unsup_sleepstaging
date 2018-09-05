%% Configuration
homedir=pwd;

SBJ_ID = 'ME_N3';
MODE = "MAIN";
%%
%NUM_OF_MAIN_CLUSTERS=2;

for NUM_OF_MAIN_CLUSTERS=2:6

% Comment the following two loops for single main clusters

if strcmp(MODE, "MAIN")
    MAX_TARGET_MAIN_CLUSTER=1;
    MAX_SUB_CLUSTERS = 2
else
    MAX_TARGET_MAIN_CLUSTER=NUM_OF_MAIN_CLUSTERS;
    MAX_SUB_CLUSTERS = 7;
end

for TARGET_MAIN_CLUSTER=1:MAX_TARGET_MAIN_CLUSTER
for NUM_OF_SUB_CLUSTERS=2:MAX_SUB_CLUSTERS
    
%%
kmeans_clustering_configuration;

%TARGET_MAIN_CLUSTER=1;
% NUM_OF_SUB_CLUSTERS=2;

if strcmp(MODE, "MAIN")
    fig_filename = sprintf('Main_%d_Clusters', NUM_OF_MAIN_CLUSTERS);
else
    fig_filename = sprintf('TotalMain_%d_Cluster_%d_SubCluster_%d', ...
        NUM_OF_MAIN_CLUSTERS, TARGET_MAIN_CLUSTER, NUM_OF_SUB_CLUSTERS);
end

PARTICIPANT=SBJ_ID;
PARTICIPANT_SECONDARY_ID='_bipolar';
FIGURE_VISIBLE = 'off';

BASE_PATH=strcat('/Volumes/Spaceship/Voss_Lucid/', PARTICIPANT, PARTICIPANT_SECONDARY_ID, '/1_CHANS');
%STAGE_LOAD_FILENAME=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, '_1_EEG_Main.mat');
% STAGE_LOAD_FILENAME=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, '_1_EEG_Main_6_Clusters.mat');

if (strcmp(MODE, "MAIN"))
    MAIN_CLUSTER_FILE_PREFIX=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, '_1_EEG_Main_', ...
        num2str(NUM_OF_MAIN_CLUSTERS), '_Clusters');
    STAGE_LOAD_FILENAME=strcat(MAIN_CLUSTER_FILE_PREFIX,  '.mat');
    STAGE_LOAD_CSV=strcat(MAIN_CLUSTER_FILE_PREFIX,  '.csv');
else
% STAGE_LOAD_FILENAME=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, '_Cluster_',num2str(TARGET_MAIN_CLUSTER), '_1_EEG_', num2str(SUB_CLUSTERS_NUM), '_REM_substages.mat');
% SUB_CLUSTER_FILE_PREFIX=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, '_Cluster_',num2str(TARGET_MAIN_CLUSTER), ...
%     '_1_EEG_', num2str(NUM_OF_SUB_CLUSTERS), '_REM_substages');
    SUB_CLUSTER_FILE_PREFIX=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, ...
        '_TotalMain_', num2str(NUM_OF_MAIN_CLUSTERS), ...
        '_Cluster_',num2str(TARGET_MAIN_CLUSTER), '_1_EEG_', num2str(NUM_OF_SUB_CLUSTERS), ...
        '_substages');
    STAGE_LOAD_FILENAME=strcat(SUB_CLUSTER_FILE_PREFIX, '.mat');
    STAGE_LOAD_CSV=strcat(SUB_CLUSTER_FILE_PREFIX, '.csv');
end

hctsafile=STAGE_LOAD_FILENAME;

PLOT_SMOOTH=0;

EDF_FILE=strcat('/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/', PARTICIPANT, '_data.mat');

%potential_lucid_dreaming_epoch_start = [17550, 20940, 21210, 22200, 22440, 22530];
potential_lucid_dreaming_epoch_start=[];
lucid_dreaming_epochs=(potential_lucid_dreaming_epoch_start./30)+1;
%%
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

channel_info=[(1:length(header.channelname))',string(header.channelname)];

kmeans_clustering_configuration_channels;

%% Pre-processing

% Wo = 50/(fs/2);  BW = Wo/35;
% [b,a] = iirnotch(Wo,BW);
% data(chans.EEG,:)=filter(b,a,data(chans.EEG,:));
% 
% % Bandpass between [0.5 70]
% [b,a]=butter(2,[0.5,70]/(fs/2),'bandpass');
% sub_ts=data(chans.EEG,:)';
% sub_ts=filtfilt(b,a,sub_ts);
% data(chans.EEG,:) = sub_ts';
% 
% [b,a]=butter(2,[0.5,15]/(fs/2),'bandpass');
% sub_ts=data(chans.EOG,:)';
% sub_ts=filtfilt(b,a,sub_ts);
% data(chans.EOG,:) = sub_ts';

%%
t=load(hctsafile);
ts1=struct2table(t.TimeSeries);
orig_datamat = t.TS_DataMat;

ts1_labels = table2array(ts1(:,1))';
ts1_algo_groups = str2num(char(table2array(ts1(:,2))));
unique_groups1=unique(ts1_algo_groups);  

% [faxis, pow] = get_PowerSpec(data(chans.EEG(1), :), 200, 1, 1); 
%figure;
set(0,'DefaultFigureVisible', FIGURE_VISIBLE);

colors=num2cell(jet(length(unique_groups1)* 5), 2);
legends = [];
legend_labels = {};
figure;

sub_clusters=[];

for i = 1:length(unique_groups1)
    
    indices = find(ts1_algo_groups==i);
    fprintf("Cluster %d has %d epochs\n", i, size(indices, 1));
    
    sub_ts1 = ts1(indices, :);
    
    legend_labels{i} = sprintf('Cluster %d (%d epochs)', i, size(sub_ts1, 1));
    
    f=[];
    p=[];
    p_ld = [];
    non_lp_count = 0;
    lp_count = 0;
    
    for j = 1:size(sub_ts1, 1)
        s = split(table2array(sub_ts1(j, 1)), '_');
        epoch_no = str2double(string(s(3, 1)));
        
        sub_clusters = [sub_clusters; i, epoch_no, (epoch_no-1)*epoch_seconds];
        
        startIndex = ((epoch_no - 1) * epoch_seconds * timeseries_sampling_rate) + 1;
        endIndex = (epoch_no * epoch_seconds * timeseries_sampling_rate);
        
        for cc = 1:length(chans.EEG)
            
            [faxis, pow] = get_PowerSpec(data(chans.EEG(cc), (startIndex:endIndex)), 200, 0, 0); 

%             faxis = faxis(faxis >= 0.5 & faxis <= 70);
%             pow = pow(find(faxis >= 0.5 & faxis <= 70));
            
            f = [f; faxis];
            
            non_lp_count = non_lp_count+1;
            p = [p; 100.*(pow./sum(pow))];
        end
    end
    
    f = mean(f);
    p = mean(p);
    
    %f=f(f<=48);
%     p(f>48) = 0;
%     p=(100*p)/sum(p);
    
    %%
    %p(faxis>48)=0;
    %p=10*log10(p);
    if (PLOT_SMOOTH == 0)
        h = plot(f, smooth(p*1000), 'Color', colors{i*4}, 'LineWidth', 2);
    else
        h = plot(f, smooth(p*1000, PLOT_SMOOTH), 'Color', colors{i*4}, 'LineWidth', 2);
    end
    
%     legends = [legends h];
%     legend_labels{i} = sprintf('Cluster %d\n', i);
%     xlabel('Frequency [Hz]')
%     ylabel('Power [dB]')
%     title(sprintf('Power Spectrum for Cluster %d\n', i));
    hold on;
end

%% Output CSV subclusters
if size(sub_clusters, 1) > 0
    sub_clusters_table = array2table(sub_clusters, 'VariableNames', {'Cluster', 'Epoch_No', 'Epoch_Start_Seconds'});
    writetable(sub_clusters_table, STAGE_LOAD_CSV);
end

set(gca, 'FontSize', 18);
xlabel('Frequency (Hz)')
ylabel('Relative Power (%)')
xlim([0.5 48]);
xticks([0:8:48]);
set(gca, 'YScale', 'log');
%ylim([0 60]);
title(sprintf('Power Spectrum (Channels C3-C4)'));
ax = gca;
%ax.TitleFontSizeMultiplier = 2;

legend(legends, legend_labels, 'FontSize', 18);
grid on;
    
hold off;

% saveas(gc f, strcat(BASE_PATH, filesep, 'Main_Cluster_', num2str(NUM_OF_MAIN_CLUSTERS), '.png'));
saveas(gcf, strcat(BASE_PATH, filesep, 'images', filesep, fig_filename), 'png');
%export_fig(strcat(BASE_PATH, filesep, 'images', filesep, fig_filename), '-eps');

end
end

%%
clearvars -except SBJ_ID MODE;

end


%%
function [faxis, pow] = get_PowerSpec(signal, SamplingRate, DecibelsFlag ,plotFlag)
% For example for a 10sec segment sampled at millisec resolution use:
% get_PowerSpec(signal, SamplingRate, DecibelsFlag ,plotFlag)
if (nargin < 4), plotFlag = 0; end
if (nargin < 3), DecibelsFlag = 0; end
% get power
pow = (abs(fft(signal)).^2)/length(signal);
% convert to decibels
if DecibelsFlag==1
    pow = 10*log10(pow/max(pow));
end
% first half of data without negative frequencies
pow = pow(1:min(floor(length(signal)/2)+1,length(signal)));
% define df and fNQ
df = 1/(length(signal)/SamplingRate);
fNQ = SamplingRate/2;

faxis = (0:df:fNQ);

if (plotFlag)
    figure;
    plot(faxis, pow);
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
end
end
