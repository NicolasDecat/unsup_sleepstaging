configuration_settings;

%% Initialise everything
homedir = pwd;
cd(HCTSA_DIR)
startup
cd(homedir)

annotation = load(ANSWER_FILE);
label = annotation.sleepstage;

fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);

feat_name = features{1,2};
% All operation names
hctsafile = HCTSA_FILE;
all_op = load(hctsafile,'Operations');

% Check operation name, get feat_id
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

%% Start of real stuff
t=load('180001_HCTSA/HCTSA_N.mat');
datamat = t.TS_DataMat;
quality = t.TS_Quality;
ts = t.TimeSeries;

single_channel_size = size(datamat,1)/3;

td=struct2table(ts);
td(:,2)=cellstr([string(label); string(label); string(label)]);
ts=table2struct(td);

% startIndex=1; % No trim
startIndex=335;

% 1 channel
t.TS_DataMat = [datamat(startIndex:single_channel_size, feat_id)];
t.TS_Quality = [quality(startIndex:single_channel_size, feat_id)];
t.TimeSeries = [ts(startIndex:single_channel_size)];

% 2 channels
% t.TS_DataMat = [datamat(startIndex:single_channel_size, feat_id) datamat(startIndex+single_channel_size:single_channel_size*2, feat_id)];
% t.TS_Quality = [quality(startIndex:single_channel_size, feat_id) quality(startIndex+single_channel_size:single_channel_size*2, feat_id)];
% t.TimeSeries = [ts(startIndex:single_channel_size);ts(startIndex+single_channel_size:single_channel_size*2)];

% 3 channels
% t.TS_DataMat = [datamat(startIndex:single_channel_size, feat_id) datamat(single_channel_size+startIndex:single_channel_size*2, feat_id) datamat(single_channel_size*2+startIndex:single_channel_size*3, feat_id)];
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

t.Operations = [reduced_operations reduced_operations reduced_operations];

save('180001_HCTSA/HCTSA_N_EEG/HCTSA.mat', '-struct', 't');
