% %% CONFIGURATION

% Global
n_cluster_samples=5;
n_clust = 16;
n_channels = 1;

kmeans_clustering_configuration;

%% Load timeseries information
ts = load(hctsafile,'TimeSeries');
ts = ts.TimeSeries;
tst = struct2table(ts);
selected_timeseries = tst(C4_COL, :);
EOG_timeseries = tst(EOG_COL,:);
EMG_timeseries = tst(EMG_COL, :);

n_epoch_ids=input(['Please enter the epoch id: [1-' num2str(size(tst, 1)/no_of_channels) ']:']);

set(gcf,'Visible','off')

for i=1:length(n_epoch_ids)
    n_epoch_id = n_epoch_ids(i);
    
    single_signal=cell2mat(table2array(selected_timeseries(n_epoch_id, 4)));
    signals{1}.signal=single_signal;
    signals{1}.signal_type='EEG';
    signals{1}.min = -1700;
    signals{1}.max = 1700;
    single_signal_length=length(single_signal);

    single_signal=cell2mat(table2array(EOG_timeseries(n_epoch_id, 4)));
    signals{2}.signal=single_signal;
    signals{2}.signal_type='EOG';
    signals{2}.min = -3950;
    signals{2}.max = 3950;

    single_signal=cell2mat(table2array(EMG_timeseries(n_epoch_id, 4)));
    signals{3}.signal=single_signal;
    signals{3}.signal_type='EMG';
    signals{3}.min = 600;
    signals{3}.max = -600;

    tslabel = string(table2cell(selected_timeseries(n_epoch_id, 2)));
    labels = split(tslabel, ',');
    timeseg = labels(2);
    timeseg = strrep(timeseg, 'timeseg_', '');
    tseg = str2double(timeseg);
    timeseg_seconds_start = (tseg-1)*epoch_seconds;    
    timeseg_seconds_end = tseg*epoch_seconds;

    draw_figure(epoch_seconds, signals, timeseries_sampling_rate, epoch_seconds, (size(ts,1)/no_of_channels)*single_signal_length, ...
        sprintf('Cluster: %d TS: %s Time: between %d secs and %d secs' , i, timeseg, timeseg_seconds_start, timeseg_seconds_end));
end

set(gcf,'Visible','on');

%             imagename = strcat('LD_Clust_',num2str(i,'%01d'),'_',num2str(j),'_',num2str(tseg,'%04d'),'.png');
%             saveas(gcf,[output_folder filesep imagename]) % saveas, imwrite or imsave? print(imagename,'-dpng')?
%             close



function draw_figure(tmax, selectedSignal, samplingRate, epoch_seconds, whole_ts_length, title_array)
    % Get number of signals
    num_signals = length(selectedSignal);
    timeID=1;
    x0=0;
    y0=0;
    width=1800;
    height=900;
    
    figure;
    for s = 1:num_signals
        % Get signal
        signal =  selectedSignal{s}.signal;
%         signal = downsample(signal, 2);
%         samplingRate = samplingRate/2;
        
        t = [0:whole_ts_length-1]/samplingRate; % = record_duration

        % Parameters for normalisation - use global max and min if amplitude
        % matters. If not, set an arbitary value
%         sigMin = -0.3; %min(signal);
%         sigMax = 0.3; %max(signal);
        sigMin = selectedSignal{s}.min;
        sigMax = selectedSignal{s}.max;
        signalRange = sigMax - sigMin;
    %     
        % Identify indexes of 30 seconds of signal according to tstart, tend
        % Otherwise, indexes = find(t<=tmax);
        tStart = find(t==(timeID-1)*epoch_seconds);
        tEnd = find(t==timeID*epoch_seconds)-1;
        indexes = tStart:1:tEnd;
        signal = signal(indexes);
        time = t(1:length(indexes)); % time = t(indexes); % Hide real time, always display 0 -30 seconds 

    %% Switch-case for num_signals
         if signalRange~= 0
            signal = signal/(sigMax-sigMin);
         end
         
    switch (num_signals)
        case 1
            % Centred around 0
            signal = signal - mean(signal);
        case {2,3}
            % Add signal below the previous one
            signal = signal - mean(signal) + (num_signals - s + 1);
            %     signal = signal + (num_signals - s + 1); % Without zero-centred
            %     signal = signal - 0.5*mean(signal) + (num_signals - s + 1);
            % Plot line dividing signals
            plot(time,s-0.5*ones(1,length(time)),'color',[0.5,0.5,0.5])
    end

    % Color code signal type - can be customised + will depends on screen setting
    switch (selectedSignal{s}.signal_type)
        case 'EEG'
            ccode = [0.1,0.5,0.8];
        case 'EOG'
            ccode = [0.1,0.5,0.3];
        case 'EMG'
            ccode = [0.8,0.5,0.2];
        case 'ECG'
            ccode = [0.8,0.1,0.2];
        otherwise
            ccode = [0.2,0.2,0.2];
    end
    set(gcf,'units','pixels','position',[x0,y0,width,height]);
    plot(time, signal,'Color',ccode);
    hold on;
end
    
    %% Plot configuration
    grid on
    ax = gca;
    fig = gcf;
    
    switch (num_signals)
        case 1
            % Set axes limits
            v = axis();
            v(1:2) = [0,tmax];
            v(3:4) = [sigMin, sigMax];
            
            axis(v);
            % Set x-axis 
            xlabel('Time(sec)')
            ax.XTick = [0:epoch_seconds];
            ax.FontSize = 10;
            % Set y-axis
            ylabel('Amplitude(\muV)')
            ax.YTick = linspace(sigMin,sigMax,epoch_seconds);

        case {2,3}
            % Set axis limits
            v = axis();
            v(1:2) = [0,tmax];
            v(3:4) = [0.5 num_signals+0.5];
            axis(v);
            % Set x axis
            xlabel('Time(sec)');
            ax = gca;
            ax.XTick = [0:epoch_seconds];
            ax.FontSize = 10;
            % Set y axis labels
            ylabel('Amplitude (mV)')
            %% Without scale
            signalLabels = cell(1,num_signals); %Revert the order such that first channel stays on top
            % for s = 1:num_signals
            %     signalLabels{num_signals-s+1} = selectedHeader(s).signal_labels;
            % end
            %ax.YTick = 1:num_signals;
            %ax.YTickLabels = signalLabels;
            %% With scale
            ax.YTick = [0.55,1,1.45,1.55,2,2.45,2.55,3,3.45];
            ax.YTickLabels = {'-300',selectedSignal{3}.signal_type,'+300','-300',selectedSignal{2}.signal_type,'+300','-300',selectedSignal{1}.signal_type,'+300'};
            ax.FontSize = 15;
            % 
            % Set figure size
            fig.Units = 'pixels';
            fig.Position = [x0,y0,width, height];
            fig.Color = [0.95 0.95 0.95];
            title(title_array);

    end

    % Reduce white space ** Can be adjusted
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + 1.1*ti(1);
    bottom = outerpos(2) + 1.1*ti(2);
    ax_width = outerpos(3) - 1.1*ti(1) - 3*ti(3);
    ax_height = outerpos(4) - 1.1*ti(2) - 4*ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    hold off;
end



