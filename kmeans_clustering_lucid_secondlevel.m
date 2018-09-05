clear feat_id feat_id_test;

% %% CONFIGURATION

% Global
n_cluster_samples=1;
substage_samples=1;
n_clust = 5;
n_channels = 1;

all_3_eegs=0;

kmeans_clustering_configuration;

SUB_CLUSTERS_NUM=2;
STAGE_LOAD_FILENAME='/Volumes/Spaceship/Voss_Lucid/KJ_N1/ALL_EEG/HCTSA_N_KJ_N1_1_EEG_7_REM_substages.mat';
SUBSTAGE_SAVE_FILENAME='HCTSA_N_KJ_N1_1_EEG_7_REM_second_level_substages.mat';
TARGET_CLUSTER_TO_SPLIT=5;

%%
all_op = load(STAGE_LOAD_FILENAME,'Operations');

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
    match_input_index = 0;
    match_test_index = 0;
    
    for i = 1:length(all_op.Operations)
        name = all_op.Operations(i).Name;
        if strcmp(op_name,name)
            match_input_index = i;
            break;
        end
    end
    
    if (match_input_index > 0)
        nn=nn+1;
        feat_id(nn) = match_input_index;
        feat(nn).id = match_input_index; % all_op.Operations(i).ID % Actual operation ID
        feat(nn).name = op_name;
    end

end
clear i n nn op_name name

%% Use feat_id to select data from full op
datamat = load(STAGE_LOAD_FILENAME,'TS_DataMat');
orig_datamat = datamat.TS_DataMat;

t = load(STAGE_LOAD_FILENAME, 'TimeSeries');
ts = struct2table(t.TimeSeries);

datamat = orig_datamat;

if n_channels == 1
    datamat = datamat;
elseif n_channels == 2
    datamat = [datamat orig_datamat(EOG_COL,feat_id)];
elseif n_channels == 3
    datamat = [datamat orig_datamat(EOG_COL,feat_id) orig_datamat(EMG_COL,feat_id)];        
end


%% Find the target cluster epochs

    idx=str2double(cellstr(table2array(ts(:, 2))));
    mkdir(output_folder);
  
    %% Lucid dreaming substage
    % Supposely the most count of potential lucid dreaming is the REM
    % stage. Use that and split the REM stage into two clusters.
    rem_stage = TARGET_CLUSTER_TO_SPLIT;
    idx_rem_stage = find(idx==rem_stage);
    
    n_clust = SUB_CLUSTERS_NUM;
    datamat=datamat(idx_rem_stage,:);
    [idx, c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);    
    
    %% Find the max cluster that has the highest silhouette values
%     max_cluster = 0;
%     max_s = 0;
%     sscore = [];
%     for i = 1:20
%         [tidx,tc, tsse] = kmeans(datamat,i,'Distance','sqeuclidean',...
%                             'Display','off','Replicates',50,'MaxIter',500);
%         s=silhouette(datamat, tidx);
%         sscore = [sscore, mean(s)];
%         if (mean(s) > max_s)
%             max_s = mean(s);
%             max_cluster = i;
%         end
%     end
% 
%     disp(['Max cluster is ' num2str(max_cluster) ' with ' num2str(max_s)]);
                    
    %% Save HTCSA_N trimmed (substages)
    t=load(STAGE_LOAD_FILENAME);
    dm = t.TS_DataMat;
    quality = t.TS_Quality;
    ts = t.TimeSeries;
    original_ts = struct2table(t.TimeSeries);

    single_channel_size = size(dm,1)/6;

    startIndex=1; % No trim

    % 1 EEG channel
    t.TS_DataMat = [dm(idx_rem_stage, feat_id)];
    t.TS_Quality = [quality(idx_rem_stage, feat_id)];
    t.TimeSeries = [ts(idx_rem_stage)];

    ts = t.TimeSeries;
    td=struct2table(ts);
    td(:,2)=cellstr([string(idx)]);
    ts=table2struct(td);
    
    t.TimeSeries = ts;
    

    %%

    % 3 EEG channels
%     t.TS_DataMat = [dm(idx_rem_stage, feat_id) dm(idx_rem_stage+single_channel_size, feat_id) dm(idx_rem_stage+(single_channel_size*2), feat_id)];
%     t.TS_Quality = [quality(idx_rem_stage, feat_id) quality(idx_rem_stage+single_channel_size, feat_id) quality(idx_rem_stage+(single_channel_size*2), feat_id)];
%     t.TimeSeries = [ts(idx_rem_stage);ts(idx_rem_stage+single_channel_size);ts(idx_rem_stage+(single_channel_size*2))];
% 
%     ts = t.TimeSeries;
%     td=struct2table(ts);
%     td(:,2)=cellstr([string(idx); string(idx); string(idx)]);
%     ts=table2struct(td);
%     
%     t.TimeSeries = ts;

    % 2 channels
    % t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(startIndex+single_channel_size:single_channel_size*2, feat_id)];
    % t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(startIndex+single_channel_size:single_channel_size*2, feat_id)];
    % t.TimeSeries = [ts(startIndex:single_channel_size);ts(startIndex+single_channel_size:single_channel_size*2)];

    % 3 channels
    % t.TS_DataMat = [dm(startIndex:single_channel_size, feat_id) dm(single_channel_size+startIndex:single_channel_size*2, feat_id) dm(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
    % t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(single_channel_size+startIndex:single_channel_size*2, feat_id) quality(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
    % t.TimeSeries = [ts(startIndex:single_channel_size);ts(single_channel_size+startIndex:single_channel_size*2);ts(single_channel_size*2+startIndex:single_channel_size*3)];

    ops=t.Operations;
    reduced_operations = [];

    feattable=struct2table(feat);
    for i = 1:size(feattable(:,1), 1)
        f_name = table2array(feattable(i,2));
        for j = 1:size(ops,1)
            op = ops(j);
            if string(op.Name) == string(cell2mat(f_name))
               reduced_operations = [reduced_operations; op]; 
            end
        end
    end

%     t.Operations = [reduced_operations reduced_operations reduced_operations];
    t.Operations = [reduced_operations];
    save(strcat('/Volumes/Spaceship/Voss_Lucid/KJ_N1/ALL_EEG/', SUBSTAGE_SAVE_FILENAME), '-struct', 't');

    
%%     
%     if (length(unique(idx(592:598))) == 1 & length(unique(idx(716:735))) ==  1 & ...
%             idx(591)~=idx(592) & idx(598) ~= idx(599) & ...
%             idx(715)~=idx(716) & idx(735) ~= idx(736))
%         disp(['Cluster ' num2str(i) ' has same consecutive clusters for lucid']);
%     end
% end

% [s,~]=silhouette(datamat, idx);
% disp(mean(s));
% n_clust=max_cluster;

%% Hierarchical clustering
% eucD=pdist(datamat, 'cosine');
% clustTreeEuc=linkage(eucD, 'average');
% 
% disp(cophenet(clustTreeEuc, eucD));
% 
% [h,nodes] = dendrogram(clustTreeEuc,0);
% h_gca = gca;
% h_gca.TickDir = 'out';
% h_gca.TickLength = [.002 0];
% h_gca.XTickLabel = [];
% 
% [h,nodes] = dendrogram(clustTreeEuc,20);
% hidx = cluster(clustTreeEuc,'criterion','distance','cutoff',.185);
% idx=hidx;
% n_clust = length(unique(hidx));

%% Plot

% markersize=12;         
% feature1 = 1;
% feature2 = 2;
% 
set(gcf,'Visible','off');

% Generate random colors
% colors = distinguishable_colors(n_clust);
% figure;
% 
% for i=1:n_clust
%     plot(datamat(idx==i,feature1),datamat(idx==i,feature2),'.','color', colors(i,:),'MarkerSize',markersize);
%     legend_str{i} = ['Cluster ' num2str(i) '(SSE: ' num2str(sse(i)) ')'];
%     hold on;
% end
% 
% dcm_obj1 = datacursormode(gcf);
% set(dcm_obj1,'UpdateFcn',{@cluster_plot_data_cursor,datamat(:, feature1:feature2),idx});
% 
% title(sprintf("KMeans Clustering - Lucid - Chan: %s", "C4"));
% plot(c(:,feature1),c(:,feature2),'kx','MarkerSize',markersize,'LineWidth',2);
% %plot(c(:,1),c(:,2),'ko','MarkerSize',12,'LineWidth',2)
% 
% legend1 = legend(string(legend_str), 'Location','NW');
% set(legend1,'Location','northeastoutside');   
% 
% set(gcf,'Visible','off')
% 
% savefig([output_folder filesep 'sp_feature1_feature2.fig']);

%% Load timeseries information
ts = load(hctsafile,'TimeSeries');
ts = ts.TimeSeries;
tst = struct2table(ts);
selected_timeseries = tst(C4_COL, :);
eeg_mat = cell2mat(table2array(selected_timeseries(:,4)));

EOG_timeseries = tst(EOG_COL,:);
eog_mat = cell2mat(table2array(selected_timeseries(:,4)));

EMG_timeseries = tst(EMG_COL, :);
emg_mat = cell2mat(table2array(selected_timeseries(:,4)));

segment_idx=[];
segment_cluster=[];
for i = 1:1
    x = length(idx);
    segment_idx = [segment_idx, [(i-1)*x+1:i*x]];
    segment_cluster = [segment_cluster, idx'];
end

segment_cluster_summary=array2table([segment_idx',(segment_idx')*epoch_seconds,segment_cluster']);
writetable(segment_cluster_summary, [output_folder filesep 'subcluster_segment_information.csv']);

%% Perform validation test (using another dataset).

set(gcf,'Visible','on')

% keywords_col = [keywords_col, array2table(idx)];
% timescale=(mod(table2array(keywords_col(:,2)), 120)*5);
% timescale(timescale==0) = 600;
% keywords_col = [keywords_col, array2table(timescale)];

% draw_figure(30, selectedSignal, 200);

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

end

function txt = cluster_plot_data_cursor(~,event_obj,myarray, indexes)
    pos = get(event_obj,'Position');
    
    idx = find((abs(myarray(:,1)-pos(1))<0.0001) & (abs(myarray(:,2)-pos(2))<0.0001));
    timeseconds = (idx-1)*30;
    txt = {['X: ',num2str(pos(1))],...
           ['Y: ',num2str(pos(2))],...
           ['Index: ',num2str(idx)],...
           ['Cluster: ', num2str(indexes(idx))],...
           ['Time: ' , num2str(timeseconds)]};
end



