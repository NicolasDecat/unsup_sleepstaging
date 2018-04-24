%% Post HCTSA v.1(After feature extraction)
% Following HCTSA analysis, prepare data for further analysis eg. cross validation 
% Post HCTSA v.1 - select features before reshaping i.e. select one feature
% = select that feature of all channels

configuration_settings

%% Normalise HCTSA.m
% Normalise to remove features containing special values (NaN, complex num)
homedir = pwd;
% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)
% Load data and normalise
cd(HCTSA_DATA_DIR)
% Normalise HCTSA.mat (remove special values)
TS_normalize('scaledRobustSigmoid',[0.8,1.0]);

% Load normalised data matrix
hctsafile = load('HCTSA_N.mat');
cd(homedir)

%% Select features for classification
% Read text file containing operations/features
fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s'); % Read names, operations, classes
fclose(fileID);

% Wanted operation names
feat_name = features{1,2};

% All operation names
all_op = hctsafile.Operations; % After loading HCTSA_N.mat - All HCTSA.mat/HCTSA_N.mat have same data structure

% Check operation name, get feature IDs
opCount = 0;
for n = 1:length(feat_name)
    op_name = char(feat_name(n));
    for i = 1:length(all_op)
        name = all_op(i).Name;
        if strcmp(op_name,name)
            opCount = opCount+1;
            feat_id(opCount) = i;
            feat(opCount).id = i; % all_op.Operations(i).ID % Actual operation ID
            feat(opCount).name = cellstr(name);
        end
    end
end
% Load data matrix from HCTSA_N.mat 
datamat = hctsafile.TS_DataMat;
[timeseries,features]=size(datamat);

% Select features - Top ops for crossval.m
hctsa_temp = datamat(:,feat_id);

%% Rearrange time series
% Re-order TS, such that each row represents featuresxchannels at each time
% segment (treating channel as feature)
delimiter =',';
Keywords = SUB_cell2cellcell({hctsafile.TimeSeries.Keywords},delimiter); % Split into sub-cells using comma delimiter

%% Obtain time label, the second element in each cell
for t = 1:length(Keywords)
    % Find which component ID the name matched
    % Due to error when creating data matrix, after ICA the channel should
    % be called component instead, hence, channel F3 means component 1.
    timelabel(t) = Keywords{t,1}(2);
end
timeuniq = unique(timelabel);

% Check size of timelabel (hence, timedata)
% number of unique label*no. channels must equal to no. of timeseries
checking = length(timeuniq)*NUM_CHANNELS;
if checking ~= timeseries
    %?
end

datamat = hctsafile.TS_DataMat;
hctsa_temp = datamat;
%% Group same time series together
n = 0; % show state of featurexchannels
if NUM_CHANNELS == 1
    hctsa_ops = hctsa_temp;
else
    for t = 1:length(timeuniq)
        n=n+1; % Essentially n = length(timeuniq)
        timename = sprintf('timeseg_%d',t);
        id = find(ismember(timelabel,timename)); %Find timeseries ID of the time label
        
        % The id could be 0 if HCTSA_Normalise remove the row.
        if isempty(id)
            continue
        end
        
        temp = hctsa_temp(id,:);
        % Reshape to make a row of the same time with all the features
        hctsa_ops(n,:) = reshape(temp,1,[]);
        % timetime{time} = timename;       
    end
end
% Each column of timedata is the featurexchannel of time segments. Size of
% timedata is [number of time segment]x[number of features x channels]