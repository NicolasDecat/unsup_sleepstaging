%% Configuration
homedir=pwd;

%%
hctsafile = '/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_1800821_HCTSA_Analysis/HCTSA_C3-A2_N.mat';
C4_COL=1:1001;
epoch_seconds=30;
no_of_channels=1;
timeseries_sampling_rate=128;

EDF_FILE='/Volumes/Spaceship/TS_ccshs-trec-1800821_7chan.mat';
ANSWER_FILE='ccshs_1800821_annot.mat';
annotation = load(strcat('/Volumes/Spaceship/ccshs_datasets', filesep, ANSWER_FILE));
annotation_label = annotation.sleepstage;

MODE = 1; % 0 - Draw kmeans clustering, 1 - Draw Expert annotations

%%
t=load(EDF_FILE);
ts=t.timeSeriesData;
labels=t.labels;

%% Extract channel names, sampling frequency from header struct
lsplit = split(t.labels', '_');
channel=unique(lsplit(:, 1), 'stable');
fs = 128;      % Assume sampling rate are the same for all channels
chans = [1];
    

%%
t=load(hctsafile);
ts1=struct2table(t.TimeSeries);
orig_datamat = t.TS_DataMat;

ts1_labels = table2array(ts1(:,1))';

if (MODE == 0)
    ts1_algo_groups = str2num(char(table2array(ts1(:,2))));
elseif (MODE == 1)
    % Need to trim based on the start and end epochs
    startIndexSplit = split(string(table2array(ts1(1, 1))), '_');
    startIndex = str2double(startIndexSplit(3, 1));
    endIndexSplit = split(string(table2array(ts1(end, 1))), '_');
    endIndex = str2double(endIndexSplit(3, 1));
    
    ts1_algo_groups = annotation_label(startIndex:endIndex, :);
end



unique_groups1=unique(ts1_algo_groups);  

%%
% [faxis, pow] = get_PowerSpec(data(chans.EEG(1), :), 200, 1, 1); 
%figure;
colors=num2cell(jet(length(unique_groups1)* 5), 2);
legends = [];
legend_labels = {};
figure;

for i = 1:length(unique_groups1)
    cluster_index = unique_groups1(i);
    indices = find(ts1_algo_groups==cluster_index);
    sub_ts1 = ts1(indices, :);
    legend_labels{i} = strcat('Cluster ', num2str(i));
    
    f=[];
    p=[];
    p_ld = [];
    non_lp_count = 0;
    lp_count = 0;
    
    for j = 1:size(sub_ts1, 1)
        s = split(table2array(sub_ts1(j, 1)), '_');
        epoch_no = str2double(string(s(3, 1)));
        
        startIndex = ((epoch_no - 1) * 30 * 200) + 1;
        endIndex = (epoch_no * 30 * 200);
        
        for cc = 1:length(chans)
            label_name=string(table2array(sub_ts1(j, 1)));
            label_idx = find(string(labels)==label_name);
            [faxis, pow] = get_PowerSpec(ts(label_idx, :), fs, 0, 0); 

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
    h = plot(f, smooth(p*1000), 'Color', colors{i*4}, 'LineWidth', 2)
    
%     legends = [legends h];
%     legend_labels{i} = sprintf('Cluster %d\n', i);
%     xlabel('Frequency [Hz]')
%     ylabel('Power [dB]')
%     title(sprintf('Power Spectrum for Cluster %d\n', i));
    hold on;
end


xlabel('Frequency (Hz)')
ylabel('Relative Power (%)')
xlim([0.5 48]);
xticks([0:8:48]);
set(gca, 'YScale', 'log');
%ylim([0 60]);
title(sprintf('Power Spectrum (Channels C3 and C4)'));
ax = gca;
ax.TitleFontSizeMultiplier = 2;

legend(legends, legend_labels);
grid on;
    
hold off;



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