function plot_cluster_fft(header, data, hctsafile, epoch_seconds, timeseries_sampling_rate, epoch_indices, epoch_labels, title_suffix)

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
    ts1=ts1(epoch_indices, :);
    
    unique_groups1=unique(epoch_labels);  

    % [faxis, pow] = get_PowerSpec(data(chans.EEG(1), :), 200, 1, 1); 
    %figure;
    colors=num2cell(jet(length(unique_groups1)* 5), 2);
    legends = [];
    legend_labels = {};
    figure;

    for i = 1:length(unique_groups1)

        sub_ts1 = ts1(find(epoch_labels==i), :);

        legend_labels{i} = strcat('Cluster ', num2str(i));

        f=[];
        p=[];
        p_ld = [];
        non_lp_count = 0;
        lp_count = 0;

        for j = 1:size(sub_ts1, 1)
            s = split(table2array(sub_ts1(j, 1)), '_');
            epoch_no = str2double(string(s(3, 1)));

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
    title(sprintf('Power Spectrum (Channels C3 and C4), %s', size(epoch_indices, 1), title_suffix));
    ax = gca;
    ax.TitleFontSizeMultiplier = 2;

    legend(legends, legend_labels);
    grid on;

    hold off;
%%
end

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
