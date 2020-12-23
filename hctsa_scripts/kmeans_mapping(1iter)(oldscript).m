% Changes from original script
% - Changed configuration settings
% - Changed feat_ID (to include all f0,664 find_top_X_features (not necessary)
% - Changed SELECT_TOP_200_FEATURES definition (to include all features)
% - epochSelectFunc (instead of epochSelectFunction, in statsout)


origlabels = [];
clusterdecision = [];
Testing_accuracy = [];
TestingAcc = [];
AUC = [];
NumIteration = [];
stgAUC = [];
Dataset = [];
NumChannels = [];
col = 1;


% at some point, change dataset: change current folder first: cd('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_001')   


% Configuration
addpath '/Users/nico/Documents/GitHub/unsup_sleepstaging';
configuration_settings;

%% All operation names
hctsafile = HCTSA_FILE;
all_op = load(hctsafile,'Operations');

clear i n nn op_name name

%% Use feat_id to select data from full op
datamat = load(hctsafile,'TS_DataMat');
datamat = datamat.TS_DataMat;
feat_id = 1:size(datamat,2);  %To include all features

[timeseries,features]=size(datamat);
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
    hctsa_ops = datamat(:,feat_id);

    % run('crossvalKR.m') 
    if conf == 'BALANCED_LABELED'
        epochSelectFunc = @epochSelect;
    elseif conf == 'UNLABELED_LABELED_A'
        epochSelectFunc = @epochSelect_unbalanced;
    elseif conf == 'UNLABELED_LABELED_B'
        epochSelectFunc = @epochSelect_unbalanced_b;
    else
        epochSelectFunc = @epochSelect;
    end
    
%     ff=load('180001_feature_corr_int_1.mat');
%     ff=ff.values;
%     ff=find(ff==1);
%     SELECT_TOP_200_FEATURES=ff';
single_channel_size=size(hctsa_ops,1)/NUM_CHANNELS;
eeg_ops=hctsa_ops(1:single_channel_size,:);
eog_ops=hctsa_ops(single_channel_size+1:single_channel_size*2,:);
emg_ops=hctsa_ops(single_channel_size*2+1:single_channel_size*3,:);
% Select only the start and end
% endW = [335,381,392,376,175];
% endS = [1373,1441,1337,1530,1491];
% START = endW(WHICH_DATA);
% END = endS(WHICH_DATA); 
% 
% N=size(hctsa_ops,2);
% 
% eeg_ops=eeg_ops(START:END, :);
% eog_ops=eog_ops(START:END, :);
% emg_ops=emg_ops(START:END, :);

% Generate all combinations of three channels
m = [];
max=200;

max_test_performance = 0;
max_combo = [];

for emg = 10:10:max
    for eog = max-emg:-10:10
        for eeg = max-emg-eog:-10:10
            e = [eeg, eog, emg];
            if sum(e) == 200
                m = [m; e];
            end
        end
    end
end

% best_perf=0;
% best_m = [];
% all_m = [];
% 
% for i=1:size(m,1)

curr_m = [60 110 30];

% eeg_values=find_top_X_features(eeg_ops, curr_m(1));
% eog_values=find_top_X_features(eeg_ops, curr_m(2))+N;
% emg_values=find_top_X_features(eeg_ops, curr_m(3))+(N*2);
% % 
% all_values=find_top_X_features([eeg_ops eog_ops emg_ops], 200);
% t=sprintf('Performance of algorithm using top %d (EEG) %d (EOG) %d (EMG) from lowest pairwise-correlation features', length(eeg_values), ...
%     length(eog_values), length(emg_values));
    
%SELECT_TOP_200_FEATURES=[eeg_values eog_values emg_values];
% SELECT_TOP_200_FEATURES=all_values;
SELECT_TOP_200_FEATURES=size(hctsa_ops,2);

    [statsOut testMat scoredTest predictTest] = cross_validation_selectivefeatures(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, SELECT_TOP_200_FEATURES);
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
    num_of_features(:) = unique(statsOut.totalFeatures);
    num_of_channels=zeros(size(statsOut.scoredTrain, 1), 1);
    num_of_channels(:) = c;
    
%     all_m = [all_m; [current_m mean(iteration_testing_accuracy)]];
%     if (mean(iteration_testing_accuracy) > best_perf)
%        best_perf = mean(iteration_testing_accuracy);
%        best_m = curr_m;
%     end
    
%     if (conf=="BALANCED_LABELED")
%         %max_test_performance = mean(iteration_testing_accuracy);
%         disp(mean(iteration_testing_accuracy));
%         plot_confusion_matrix(c, statsOut.scoredTrain, statsOut.predictTrain, ...
%             statsOut.scoredTest, statsOut.predictTest, CM_SAVE_DIR, strcat(conf, num2str(c)));
% %         plot_confusion_matrix(statistics(idx).id, statistics(idx).scoredTrain, statistics(idx).svmPredictTrain, ...
% %              statistics(idx).scoredTest, statistics(idx).svmPredictTest, CM_SAVE_DIR)
% 
%     end
    
    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Unsupervised_', conf);
    row = [types, iteration', iteration_training_accuracy, iteration_testing_accuracy, num_of_features, num_of_channels];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Supervised_', conf);
    row = [types, iteration', iteration_svm_training_accuracy, iteration_svm_testing_accuracy, num_of_features, num_of_channels];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

% end
    
end
end
end

%% Draw the confusion matrix for the repeat that has maximum trainCorrect
if PLOT_CONFUSION_MATRIX
    for idx = 1:length(exps)
        plot_confusion_matrix(statistics(idx).id, statistics(idx).scoredTrain, statistics(idx).predictTrain, ...
            statistics(idx).scoredTest, statistics(idx).predictTest, CM_SAVE_DIR)
    end
end

set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying
s=save_stats; [UA, ~, idx] = unique(s(:,[1 6]));NEW_A = [UA,array2table(accumarray(idx,double(table2array(s(:,4))),[],@mean))]; NEW_A;

%
figure
x=categorical(cellstr(char(strrep(table2array(NEW_A(:,1)), '_', ' '))));
y=table2array(NEW_A(:,3))*100;
bar(x, y);
% title(t);
title('Testing accuracy')
grid on;
ylim([0 100]);
labels = arrayfun(@(value) num2str(value,'%2.1f'), y,'UniformOutput',false);
text(x,y,labels,...
  'HorizontalAlignment','center',...
  'VerticalAlignment','bottom') 

% save(strcat(CM_SAVE_DIR, filesep, 'SUP_UNSUP_HCTSA_200_DS1.mat'), 'save_stats');

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

% Testing_accuracy = table2array(NEW_A(2,3))         

run('type1auc.m')




       