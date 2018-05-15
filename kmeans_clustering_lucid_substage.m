
%% CONFIGURATION

% Global
n_cluster_samples=10;
n_clust = 3; % for substages

kmeans_clustering_configuration;

substage_output_folder='ME_N3_5_clusters';
likely_lucid_cluster=3;
include_epoch_number_more_than=5;
cluster_csv='cluster_segment_information.csv';

%%
all_op = load(hctsafile,'Operations');
OPS_FILE='reduced_ops.txt';
%%
currentdir = pwd;

fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);

%% Wanted operation names
feat_name = features{1,2};

%% Check operation name, get feat_id
nn=0;
for n = 1:length(feat_name)
    op_name = char(feat_name(n));
    for i = 1:length(all_op.Operations)
        name = all_op.Operations(i).Name;
        if strcmp(op_name,name)
            nn=nn+1;
            feat_id(nn) = i;
            feat(nn).id = i; % all_op.Operations(i).ID % Actual operation ID
            feat(nn).name = cellstr(name);
        end
    end
end
clear i n nn op_name name

%% Use feat_id to select data from full op
datamat = load(hctsafile,'TS_DataMat');
orig_datamat = datamat.TS_DataMat;

datamat = orig_datamat(C4_COL,:);
datamat = [datamat orig_datamat(EOG_COL,:)];
datamat = datamat(:,feat_id);

%% Here we read the whole cluster CSV and extract only the likelihood cluster for substage clustering (30 seconds).
cluster_table = readtable([substage_output_folder filesep cluster_csv]);
likely_cluster = cluster_table(table2array(cluster_table(:,3))==likely_lucid_cluster,:);
lucid_datamat = [];

count = 1;
temp_datamat = datamat(1,:); % starts with first row
datamat_indexes = [];
for i = 2:size(likely_cluster,1)
   temp_datamat = [temp_datamat; datamat(i,:)];
   
   if (table2array(likely_cluster(i,1))==table2array(likely_cluster(i-1,1))+1)
      count = count + 1; 
   else
      if (count >= include_epoch_number_more_than)
        cluster_csv_index = i-count:i-1;
        
        datamat_index = table2array(likely_cluster(cluster_csv_index,1));
        datamat_indexes = [datamat_indexes datamat_index'];
        lucid_datamat = [lucid_datamat; datamat(datamat_index, :)];
      end
      
      count = 1;
      temp_datamat = datamat(i, :);
   end
end

lucid_datamat=lucid_datamat(:,1:16);

%% Perform clustering
[idx,c, sse] = kmeans(lucid_datamat,n_clust,'Distance','sqeuclidean',...
                    'Display','off','Replicates',50,'MaxIter',500);

segment_idx=[];
segment_cluster=[];
for i = 1:1
    x = length(idx);
    segment_idx = datamat_indexes;
    segment_cluster = [segment_cluster, idx'];
end

segment_cluster_summary=array2table([segment_idx',(segment_idx')*30,segment_cluster']);
writetable(segment_cluster_summary, [substage_output_folder filesep 'sub_cluster_segment_information.csv']);


%% Plot
% markersize=12;         
% 
% for i=2:16
% 
%     figure;
%     
%     feature1 = i-1;
%     feature2 = i;
%     
%     plot(datamat(idx==1,feature1),datamat(idx==1,feature2),'r.','MarkerSize',markersize)
%     hold on;
%     plot(datamat(idx==2,feature1),datamat(idx==2,feature2),'b.','MarkerSize',markersize)
%     plot(datamat(idx==3,feature1),datamat(idx==3,feature2),'k.','MarkerSize',markersize)
% %     plot(datamat(idx==4,feature1),datamat(idx==4,feature2),'g.','MarkerSize',markersize)
% %     plot(datamat(idx==5,feature1),datamat(idx==5,feature2),'c.','MarkerSize',markersize)
% %     plot(datamat(idx==6,feature1),datamat(idx==6,feature2),'m.','MarkerSize',markersize)
% 
%     title(sprintf("KMeans Clustering - Lucid - Chan: %s F: %d-%d", "C4", feature1, feature2));
%     plot(c(:,feature1),c(:,feature2),'kx','MarkerSize',markersize,'LineWidth',2)
%     %plot(c(:,1),c(:,2),'ko','MarkerSize',12,'LineWidth',2)
% 
%     legend1 = legend(['Cluster 1 (SSE: ' num2str(sse(1))],...
%            ['Cluster 2 (SSE: ' num2str(sse(2))],...
%            ['Cluster 3 (SSE: ' num2str(sse(3))],...
%            'Centroids', 'Location','NW');
% %            ['Cluster 4 (SSE: ' num2str(sse(4))],...
% %            ['Cluster 5 (SSE: ' num2str(sse(5))],...
% %            ['Cluster 6 (SSE: ' num2str(sse(6))],...
% 
%     set(legend1,'Location','northeastoutside');   
%     hold off;
% end
%% Load timeseries information

% ts = load(hctsafile,'TimeSeries');
% ts = ts.TimeSeries;
% tst = struct2table(ts);
% selected_timeseries = tst(C4_COL, :);
% eeg_mat = cell2mat(table2array(selected_timeseries(:,4)));
% 
% EOG_timeseries = tst(EOG_COL,:);
% eog_mat = cell2mat(table2array(selected_timeseries(:,4)));
% 
% EMG_timeseries = tst(EMG_COL, :);
% emg_mat = cell2mat(table2array(selected_timeseries(:,4)));
% 
% set(gcf,'Visible','off')
% 
% for i = 1:n_clust
%     selected_eeg_cluster = selected_timeseries(datamat_indexes(idx==i),:);
%     selected_eog_cluster = EOG_timeseries(datamat_indexes(idx==i),:);
%     selected_emg_cluster = EMG_timeseries(datamat_indexes(idx==i),:);
%     random_cluster_ts_idx = randperm(size(selected_eeg_cluster, 1));
%     
%     if (length(random_cluster_ts_idx) <= n_cluster_samples)
%         random_cluster_ids = random_cluster_ts_idx(1:length(random_cluster_ts_idx));
%     else
%         random_cluster_ids = random_cluster_ts_idx(1:n_cluster_samples);
%     end
%     
%     for j = 1:length(random_cluster_ids)
%         random_cluster_id = random_cluster_ids(j);
% 
%         single_signal=cell2mat(table2array(selected_eeg_cluster(random_cluster_id, 4)));
%         signals{1}.signal=single_signal;
%         signals{1}.signal_type='EEG';
%         signals{1}.min = -1700;
%         signals{1}.max = 1700;
%         single_signal_length=length(single_signal);
% 
%         single_signal=cell2mat(table2array(selected_eog_cluster(random_cluster_id, 4)));
%         signals{2}.signal=single_signal;
%         signals{2}.signal_type='EOG';
%         signals{2}.min = -3950;
%         signals{2}.max = 3950;
% 
%         single_signal=cell2mat(table2array(selected_emg_cluster(random_cluster_id, 4)));
%         signals{3}.signal=single_signal;
%         signals{3}.signal_type='EMG';
%         signals{3}.min = 600;
%         signals{3}.max = -600;
% 
%         tslabel = string(table2cell(selected_eeg_cluster(random_cluster_id, 2)));
%         labels = split(tslabel, ',');
%         timeseg = labels(2);
%         timeseg = strrep(timeseg, 'timeseg_', '');
%         tseg = str2double(timeseg);
%         timeseg_seconds_start = (tseg-1)*30;    
%         timeseg_seconds_end = tseg*30;
%         
%         set(gcf,'Visible','off')
%         draw_figure(epoch_seconds, signals, timeseries_sampling_rate, (size(ts,1)/no_of_channels)*single_signal_length, ...
%             sprintf('Cluster: %d TS: %s Time: between %d secs and %d secs' , i, timeseg, timeseg_seconds_start, timeseg_seconds_end));
%         
%         set(gcf,'Visible','off');
%         
%         imagename = strcat('LD_SubStage_Clust_',num2str(i,'%01d'),'_',num2str(j),'_',num2str(tseg,'%04d'),'.png');
%         saveas(gcf,[substage_output_folder filesep imagename]) % saveas, imwrite or imsave? print(imagename,'-dpng')?
%         close
%     end
% end

set(gcf,'Visible','on')

% keywords_col = [keywords_col, array2table(idx)];
% timescale=(mod(table2array(keywords_col(:,2)), 120)*5);
% timescale(timescale==0) = 600;
% keywords_col = [keywords_col, array2table(timescale)];

% draw_figure(30, selectedSignal, 200);

function draw_figure(tmax, selectedSignal, samplingRate, whole_ts_length, title_array)
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
        tStart = find(t==(timeID-1)*30);
        tEnd = find(t==timeID*30)-1;
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
            ax.XTick = [0:30];
            ax.FontSize = 10;
            % Set y-axis
            ylabel('Amplitude(\muV)')
            ax.YTick = linspace(sigMin,sigMax,30);

        case {2,3}
            % Set axis limits
            v = axis();
            v(1:2) = [0,tmax];
            v(3:4) = [0.5 num_signals+0.5];
            axis(v);
            % Set x axis
            xlabel('Time(sec)');
            ax = gca;
            ax.XTick = [0:30];
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

end




