% Select features from reduced_ops.txt, make another .mat file to work on
% cross-validation.
configuration_settings;

homedir = pwd;
%% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)

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
single_channel_size = timeseries/3;

% CONFIG
hctsa_ops = datamat(1:single_channel_size,feat_id);
%hctsa_ops = datamat(single_channel_size + 1 : single_channel_size*2,feat_id);
%hctsa_ops = datamat(single_channel_size*2 + 1 : single_channel_size*3,feat_id);

%% Run cross-validation code
% Change the number of operations
set(0,'DefaultFigureVisible','off') % Remove this to disable the figure displaying (sometimes it could be lots of figures!)
exps = EXPS_TO_RUN; % This allow us to selectively choose which experiment to run
statistics = [];

save_stats_columns = {'Type', 'Iteration', 'TrainingAccuracy', 'TestingAccuracy', 'FeatureId', 'FeatureName','Channels'};
save_stats = array2table(zeros(0,length(save_stats_columns)));
save_stats.Properties.VariableNames = save_stats_columns;
    
for c = NUM_CHANNELS_TO_RUN
for configuration = 1:length(CONFIGURATIONS_TO_RUN)
    conf = CONFIGURATIONS_TO_RUN(configuration);

for l  = 1:size(feat_id,2)
    k = feat_id(l);
    
    hctsa_ops = datamat(:, k);
    
    % NOTE: Here the number of channel is always 1 because we use
    % exclusively the data from this channel.
    statsOut = cross_validation(k, hctsa_ops, CM_SAVE_DIR, 1, @epochSelect);
    
    [~, statsOut.complexity]=size(hctsa_ops);
    %statsOut.complexity = k;
    statsOut.id = k;
    statistics = [statistics, statsOut];

    iteration=1:size(statsOut.scoredTrain, 1);
    iteration_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.predictTrain)'))/size(statsOut.scoredTrain, 2))';
    iteration_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.predictTest)'))/size(statsOut.scoredTest, 2))';
    iteration_svm_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.svmPredictTrain)'))/size(statsOut.scoredTrain, 2))';
    iteration_svm_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.svmPredictTest)'))/size(statsOut.scoredTest, 2))';
    num_of_features=zeros(size(statsOut.scoredTrain, 1), 1);
    num_of_features(:) = size(hctsa_ops, 2);
    num_of_channels=zeros(size(statsOut.scoredTrain, 1), 1);
    num_of_channels(:) = 1; % NOTE: Here hardcoded to 1
    
    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Unsupervised_', conf);
    
    feature_ids = ones(size(iteration, 2), 1) .* k;
    feature_names = cell(1, size(iteration, 2));
    feature_names(:) = cellstr(all_op.Operations(k).Name);

    channels = cell(1, size(iteration, 2));
    
    if (c == 1)
        channels(:) = cellstr("EEG");
    elseif(c==2)
        channels(:) = cellstr("EOG");
    elseif(c==3)
        channels(:) = cellstr("EMG");
    end

    row = [types, iteration', double(iteration_training_accuracy), double(iteration_testing_accuracy), feature_ids, string(cell2mat(feature_names')), channels'];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Supervised_', conf);
    row = [types, iteration', double(iteration_svm_training_accuracy), double(iteration_svm_testing_accuracy), feature_ids, string(cell2mat(feature_names')), channels'];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];
    
end
end
end

%% Draw the confusion matrix for the repeat that has maximum trainCorrect
if PLOT_CONFUSION_MATRIX
    for idx = 1:size(feat_id,2)
        plot_confusion_matrix(statistics(idx).id, statistics(idx).scoredTrain, statistics(idx).predictTrain, ...
            statistics(idx).scoredTest, statistics(idx).predictTest, CM_SAVE_DIR)
    end
end

set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying

save(strcat(CM_SAVE_DIR, filesep, 'SINGLE_FEATURES_DS1.mat'), 'save_stats');

%% Plot output accuracy
accuracy_train = [];
accuracy_test = [];
complexity = [];

for l=1:size(feat_id,2)
    accuracy_train = [accuracy_train, statistics(l).output.trainCorrect];
    accuracy_test = [accuracy_test, statistics(l).output.testCorrect];
    complexity = [complexity, statistics(l).complexity];
end

if PLOT_ACCURACY_REPORT
    figure;
    semilogx(complexity,accuracy_train,complexity,accuracy_test)
    legend('Training','Test')
    ylabel('Accuracy [0-1]')
    xlabel('Number of features')

    saveas(gcf, strcat(CM_SAVE_DIR, filesep, 'ACCURACY_REPORT.png'));
end