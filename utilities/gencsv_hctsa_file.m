%%
% This is a script that converts HCTSA to a CSV file (with the labels) so
% that it can be imported to WEKA tool easily.

%% Configuration
HCTSA_FILE='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/sleep_org_180001_HCTSA/HCTSA_N.mat';
ANSWER_FILE='/Volumes/Spaceship/ccshs_datasets/ccshs_1800001_annot.mat';
OPS_FILE='reduced_ops.txt';
OUTPUT_DIR='/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results/weka_datasets';
OUTPUT_FILE_PREFIX='180001_HCTSA_N';
OUTPUT_FILENAME_SUFFIX=["EEG", "EOG", "EMG"];

START=335;
END=1373;
CHANNEL_SIZE=3;

homedir = pwd;
% %% Start HCTSA tools
% cd(HCTSA_DIR)
% startup
% cd(homedir)

%% Read text file
fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);

%% Wanted operation names
feat_name = features{1,2};

%% All operation names
hctsafile = HCTSA_FILE;
all_op = load(hctsafile,'Operations');

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
datamat = datamat.TS_DataMat;

[timeseries,features]=size(datamat);
single_channel_size = timeseries/CHANNEL_SIZE;

% CONFIG
annotation = load(ANSWER_FILE);
label = annotation.sleepstage;
label=label+1;
label(label==6)=5;
label = label(START:END);
stgLab = {'W','N1','N2','N3','R'};

featt=struct2table(feat);
total_column_names = [];
hctsa_all = [];

for i = 1:CHANNEL_SIZE
    hctsa_ops = datamat((i-1)*single_channel_size+1:single_channel_size*i,feat_id);
    hctsa_ops = hctsa_ops(START:END, :);
    hctsa_all = [hctsa_all hctsa_ops];
    hctsa_ops = [hctsa_ops string((stgLab(label)'))];

    hctsa_table= array2table(hctsa_ops);
    
    column_names = cellstr([strcat(OUTPUT_FILENAME_SUFFIX(i), '_', string(table2cell(featt(:,1)))')]);
    hctsa_table.Properties.VariableNames = [column_names cellstr("labels")];
    writetable(hctsa_table, strcat(OUTPUT_DIR, filesep, OUTPUT_FILE_PREFIX, '_', OUTPUT_FILENAME_SUFFIX(i), '.csv'));
    
    total_column_names = [total_column_names column_names];
end

hctsa_all = [hctsa_all string((stgLab(label)'))];
hctsa_table= array2table(hctsa_all);

hctsa_table.Properties.VariableNames = [total_column_names cellstr("labels")];
writetable(hctsa_table, strcat(OUTPUT_DIR, filesep, OUTPUT_FILE_PREFIX, '_ALL', '.csv'));

