function plot_power_spectrum(hctsafile, header, data, epoch_seconds, ...
    timeseries_sampling_rate, chans, STAGE_LOAD_CSV, plot_smooth)
    %% Extract channel names, sampling frequency from header struct
    channel = cellstr(header.channelname);  % Channel names + converted into cell array
    %fs = header.samplerate(1);      % Assume sampling rate are the same for all channels
    %recordtime = header.records;    % aaTotal recorded time (seconds)
    %nchannels = header.channels;    % Number of channels

    % Amplitude calibration
    phys_range = header.physmax - header.physmin;
    dig_range = header.digimax - header.digimin;

    gain = phys_range ./ dig_range;

    calibrate_data=data;
    for i = 1:size(calibrate_data, 1)
       calibrate_data(i, :) = (calibrate_data(i, :) - header.digimin(i)) * gain(i) + header.physmin(i);
    end
 
    % channel_info=[(1:length(header.channelname))',string(header.channelname)];
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
    ts1_algo_groups = str2num(char(table2array(ts1(:,2))));
    unique_groups1=unique(ts1_algo_groups);  

    set(0,'DefaultFigureVisible', 'off');

    %colors=num2cell(jet(length(unique_groups1)* 5), 2);
    colors=GiveMeColors(length(unique_groups1));
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
        non_lp_count = 0;

        for j = 1:size(sub_ts1, 1)
            s = split(table2array(sub_ts1(j, 1)), '_');
            epoch_no = str2double(string(s(3, 1)));

            sub_clusters = [sub_clusters; i, epoch_no, (epoch_no-1)*epoch_seconds];

            startIndex = ((epoch_no - 1) * epoch_seconds * timeseries_sampling_rate) + 1;
            endIndex = (epoch_no * epoch_seconds * timeseries_sampling_rate);

            for cc = 1:length(chans.EEG)

                [faxis, pow] = get_PowerSpec(calibrate_data(chans.EEG(cc), (startIndex:endIndex)), 200, 0, 0); 

                f = [f; faxis];

                non_lp_count = non_lp_count+1;
                p = [p; 100.*(pow./sum(pow))];
            end
        end

        f = mean(f);
        p = mean(p);

        if (plot_smooth == 0)
            h = plot(f, smooth(p*1000), 'Color', colors{i}, 'LineWidth', 2);
        else
            h = plot(f, smooth(p*1000, PLOT_SMOOTH), 'Color', colors{i}, 'LineWidth', 2);
        end
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
    title(sprintf('Power Spectrum (Channels %s)', join(string(channel(chans.EEG))', ',')));
    ax = gca;
    %ax.TitleFontSizeMultiplier = 2;

    legend(legends, legend_labels, 'FontSize', 18);
    grid on;
    hold off;
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
