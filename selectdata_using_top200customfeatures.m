% Select features from reduced_ops.txt, make another .mat file to work on
% cross-validation.
configuration_settings;

homedir = pwd;
%% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)

%% Read the feature tables
fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);
%% Wanted operation names
feat_name = features{1,2};

%% All operation names
hctsafile = HCTSA_FILE;
all_op = load(hctsafile,'Operations');

single_feat_stats = load(SINGLE_FEATURE_RESULT, 'save_stats');
single_feat_stats = single_feat_stats.save_stats;

%% Use feat_id to select data from full op
datamat = load(hctsafile,'TS_DataMat');
datamat = datamat.TS_DataMat;

[timeseries,features]=size(datamat);
max_test_performance = 0;

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
 
hctsa_ops = datamat(:,feat_id);

%% Run cross-validation code
% Change the number of operations
set(0,'DefaultFigureVisible','off') % Remove this to disable the figure displaying (sometimes it could be lots of figures!)
exps = EXPS_TO_RUN; % This allow us to selectively choose which experiment to run
statistics = [];

save_stats_columns = {'Type', 'Iteration', 'TrainingAccuracy', 'TestingAccuracy', 'NumberOfFeatures','NumberOfChannels'};
save_stats = array2table(zeros(0,length(save_stats_columns)));
save_stats.Properties.VariableNames = save_stats_columns;
    
for c = NUM_CHANNELS_TO_RUN
for configuration = 1:length(CONFIGURATIONS_TO_RUN)
    conf = CONFIGURATIONS_TO_RUN(configuration);

for l  = 1:length(exps) 
    k = exps(l); % k is the condition to select operation

    if conf == 'BALANCED_LABELED'
        epochSelectFunc = @epochSelect;
    elseif conf == 'UNLABELED_LABELED_A'
        epochSelectFunc = @epochSelect_unbalanced;
    elseif conf == 'UNLABELED_LABELED_B'
        epochSelectFunc = @epochSelect_unbalanced_b;
    else
        epochSelectFunc = @epochSelect;
    end
    
    %%
    % Construct hctsa_ops based on channels
    
    % Use unsupervised single features metrics as the basis
    ftable = single_feat_stats(startsWith(single_feat_stats.Type, 'Unsupervised') & endsWith(single_feat_stats.Type, conf), :);
    [UA, ~, idx] = unique(ftable(:,[1 5 6 7]));
    ftable = [UA,array2table(accumarray(idx,double(table2array(ftable(:,4))),[],@mean))];
    ftable = sortrows(ftable, 'Var1', 'descend');
    
    ftable_cols=[2, 3, 4];
    
    eeg_features = ftable(ftable.Channels=="EEG", ftable_cols);
    eeg_eog_features = ftable(ftable.Channels=="EEG" | ftable.Channels=="EOG", ftable_cols);
    all_features = ftable(:, ftable_cols);

    if c == 1
        top_200_features = eeg_features(1:198, :);
    elseif c == 2
        top_200_features = eeg_eog_features(1:198, :);
    elseif c == 3
         top_200_features = all_features(1:198, :);
    end

    eeg_features = table2array(top_200_features(top_200_features.Channels=="EEG", [1, 2]));
    eog_features = table2array(top_200_features(top_200_features.Channels=="EOG", [1, 2]));
    emg_features = table2array(top_200_features(top_200_features.Channels=="EMG", [1, 2]));

    eeg_top_200_features = str2double(eeg_features(:, 1)');
    eog_top_200_features = str2double(eog_features(:, 1)');
    emg_top_200_features = str2double(emg_features(:, 1)');

    [~, eeg_top_200_feature_indexes] = intersect(feat_id, eeg_top_200_features, 'stable');
    [~, eog_top_200_feature_indexes] = intersect(feat_id, eog_top_200_features, 'stable');
    [~, emg_top_200_feature_indexes] = intersect(feat_id, emg_top_200_features, 'stable');
    
    % Run the cross-validation
    statsOut = cross_validation_customfeatures(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, eeg_top_200_feature_indexes, eog_top_200_feature_indexes, emg_top_200_feature_indexes);
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
    num_of_channels(:) = c;
    
%     if (mean(iteration_testing_accuracy) > max_test_performance && conf=="BALANCED_LABELED")
    if (conf=="BALANCED_LABELED")
        %max_test_performance = mean(iteration_testing_accuracy);
        disp(mean(iteration_testing_accuracy));
        
        plot_confusion_matrix(c, statsOut.scoredTrain, statsOut.predictTrain, ...
            statsOut.scoredTest, statsOut.predictTest, CM_SAVE_DIR, strcat(conf, num2str(c)));
%         plot_confusion_matrix(statistics(idx).id, statistics(idx).scoredTrain, statistics(idx).svmPredictTrain, ...
%             statistics(idx).scoredTest, statistics(idx).svmPredictTest, CM_SAVE_DIR)

    end
    
    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Unsupervised_', conf);
    row = [types, iteration', iteration_training_accuracy, iteration_testing_accuracy, num_of_features, num_of_channels];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Supervised_', conf);
    row = [types, iteration', iteration_svm_training_accuracy, iteration_svm_testing_accuracy, num_of_features, num_of_channels];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];
    
end
end
end

%% Draw the confusion matrix for the repeat that has maximum trainCorrect
if PLOT_CONFUSION_MATRIX
    for idx = 1:c
        plot_confusion_matrix(idx, statistics(idx).scoredTrain, statistics(idx).predictTrain, ...
            statistics(idx).scoredTest, statistics(idx).predictTest, CM_SAVE_DIR)
%         plot_confusion_matrix(statistics(idx).id, statistics(idx).scoredTrain, statistics(idx).svmPredictTrain, ...
%             statistics(idx).scoredTest, statistics(idx).svmPredictTest, CM_SAVE_DIR)
    end
end



set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying
save_stats
%save(strcat(CM_SAVE_DIR, filesep, 'SUP_UNSUP_Single_Feature_200_DS1.mat'), 'save_stats');

%% Plot output accuracy
accuracy_train = [];
accuracy_test = [];
complexity = [];

for l=1:length(exps)
    k = exps(l);
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

disp('FINAL----------');
disp(max_test_performance);
% disp(max_svm_test_performance);
