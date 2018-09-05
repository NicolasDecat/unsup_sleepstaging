clear feat_id feat_id_test;

% %% CONFIGURATION

% Global
%MAIN_CLUSTERES_NUM=3;
SUB_CLUSTERS_NUM=[2:7];
%TARGET_CLUSTER_TO_SPLIT=3;

for MAIN_CLUSTERES_NUM=2:6
    
for TARGET_CLUSTER_TO_SPLIT=1:MAIN_CLUSTERES_NUM    

n_channels = 1;
all_3_eegs=0;

kmeans_clustering_configuration;

PARTICIPANT='ME_N3';
PARTICIPANT_SECONDARY_ID='_bipolar';
BASE_PATH=strcat('/Volumes/Spaceship/Voss_Lucid/', PARTICIPANT, PARTICIPANT_SECONDARY_ID, '/1_CHANS');
EDF_FILE=strcat('/Volumes/Spaceship/Voss_Lucid/edf_init_mat_files/', PARTICIPANT, '_data.mat');
%STAGE_LOAD_FILENAME=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, '_1_EEG_Main.mat');
STAGE_LOAD_FILENAME=strcat(BASE_PATH, filesep, 'HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, ...
    '_1_EEG_Main_', num2str(MAIN_CLUSTERES_NUM), '_Clusters.mat');

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

if all_3_eegs
    O2_COL=C4_COL+length(C4_COL);
    F4_COL=O2_COL+length(C4_COL);
    datamat = [orig_datamat(C4_COL, feat_id) orig_datamat(O2_COL, feat_id) orig_datamat(F4_COL, feat_id)];
else
    datamat = orig_datamat(C4_COL, feat_id);
end

if n_channels == 1
    datamat = datamat;
elseif n_channels == 2
    datamat = [datamat orig_datamat(EOG_COL,feat_id)];
elseif n_channels == 3
    datamat = [datamat orig_datamat(EOG_COL,feat_id) orig_datamat(EMG_COL,feat_id)];        
end

    idx=str2double(cellstr(table2array(ts(:, 2))));
    mkdir(output_folder);
  
    %% Lucid dreaming substage
    % Supposely the most count of potential lucid dreaming is the REM
    % stage. Use that and split the REM stage into two clusters.
%     rem_stage = mode(potential_lucid_dreaming_stages);
    rem_stage = TARGET_CLUSTER_TO_SPLIT;
    idx_rem_stage = find(idx==rem_stage);
    datamat=datamat(idx_rem_stage,:);
    
    calculate_opt_k(datamat, "Sub-cluster - ");
    
    for sc = 1:size(SUB_CLUSTERS_NUM, 2)
        subcluster=SUB_CLUSTERS_NUM(sc);
        SUBSTAGE_SAVE_FILENAME=strcat('HCTSA_N_', PARTICIPANT, PARTICIPANT_SECONDARY_ID, ...
            '_TotalMain_', num2str(MAIN_CLUSTERES_NUM), ...
            '_Cluster_', num2str(TARGET_CLUSTER_TO_SPLIT), '_1_EEG_', num2str(subcluster), ...
            '_substages.mat');

        n_clust = subcluster;
        [idx, c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                            'Display','off','Replicates',50,'MaxIter',500);    

    %     disp(idx(find(ismember(idx_rem_stage, (592:598)))));

        %% Save HTCSA_N trimmed (substages)
        t=load(hctsafile);
        dm = t.TS_DataMat;
        quality = t.TS_Quality;
        ts = t.TimeSeries;
        original_ts = struct2table(t.TimeSeries);

        single_channel_size = size(dm,1)/no_of_channels;

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
        save(strcat(BASE_PATH, filesep, SUBSTAGE_SAVE_FILENAME), '-struct', 't');
    end
    

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

segment_idx=[];
segment_cluster=[];
for i = 1:1
    x = length(idx);
    segment_idx = [segment_idx, [(i-1)*x+1:i*x]];
    segment_cluster = [segment_cluster, idx'];
end

segment_cluster_summary=array2table([segment_idx',(segment_idx')*epoch_seconds,segment_cluster']);
writetable(segment_cluster_summary, [output_folder filesep 'subcluster_segment_information.csv']);

end
end

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



