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
    if k==0
        hctsa_ops = datamat(:,feat_id(1:1));
    elseif k==1
        hctsa_ops = datamat(:,feat_id(1:10));
    elseif k==2
        hctsa_ops = datamat(:,feat_id(1:25));
    elseif k==3
        hctsa_ops = datamat(:,feat_id(1:50));
    elseif k==4
        hctsa_ops = datamat(:,feat_id(1:100));
    elseif k==5 % Top 198 features
        disp(size(feat_id));
        hctsa_ops = datamat(:,feat_id);
        %hctsa_ops = datamat(:,feat_id(1:200));
    elseif k==6 % Random 500 features
        rand_id = randperm(features,500);
        hctsa_ops = datamat(:,rand_id);
    elseif k==7
        rand_id = randperm(features,1000);
        hctsa_ops = datamat(:,rand_id);
    elseif k==8
        rand_id = randperm(features,2000);
        hctsa_ops = datamat(:,rand_id);
    elseif k==9
        rand_id = randperm(features,5000);
        hctsa_ops = datamat(:,rand_id);
    else
        hctsa_ops = data
        
        p  
    %% Evaluate the clustering as a whole
    clust = zeros(size(datamat,1), 10);
    for i = 1:10
        clust(:, i) = kmeans(datamat, i,'Distance','sqeuclidean',...
                            'Display','off','Replicates',50,'MaxIter',500);
    end

    s=evalclusters(datamat, clust, 'silhouette');
    c=evalclusters(datamat, clust, 'CalinskiHarabasz');
    db=evalclusters(datamat, clust, 'DaviesBouldin');

    fprintf("Full dataset all features: Silhouette: %.04f (%d) CalinskiHarabasz: %.04f (%d) DaviesBouldin: %.04f (%d)\n\n", ...
        s.CriterionValues(s.OptimalK), ...
        s.OptimalK, ...
        c.CriterionValues(c.OptimalK), ...
        c.OptimalK, ...
        db.CriterionValues(db.OptimalK), ...
        db.OptimalK);
        
    
    % run('crossvalKR.m') 
    if conf == 'BALANCED_LABELED'
        epochSelectFunc = @epochSelect;
    elseif conf == 'UNBALANCED_LABELED_A'
        epochSelectFunc = @epochSelect_unbalanced;
    elseif conf == 'UNBALANCED_LABELED_B'
        epochSelectFunc = @epochSelect_unbalanced_b;
    else
        warning(strcat(conf, ' does not match any of the configuration, please check the code...'));
        epochSelectFunc = @epochSelect;
    end
    
    %%
    statsOut = cross_validation(k, hctsa_ops, CM_SAVE_DIR, c, NUM_CHANNELS, epochSelectFunc);
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
    
    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Unsupervised_', conf);
    row = [types, iteration', iteration_training_accuracy, iteration_testing_accuracy, num_of_features, num_of_channels];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

    types=strings(size(statsOut.scoredTrain, 1), 1);
    types(:)=strcat('Supervised_', conf);
    row = [types, iteration', iteration_svm_training_accuracy, iteration_svm_testing_accuracy, num_of_features, num_of_channels];
    save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];
    
    if PLOT_CONFUSION_MATRIX
        plot_confusion_matrix(strcat('_chan_', num2str(c), '_', conf), statsOut.scoredTrain, statsOut.predictTrain, ...
            statsOut.scoredTest, statsOut.predictTest, CM_SAVE_DIR)
    end
    
end
end
end

%% Draw the confusion matrix for the repeat that has maximum trainCorrect
% if PLOT_CONFUSION_MATRIX
%     for idx = 1:length(exps)
%         plot_confusion_matrix(statistics(idx).id, statistics(idx).scoredTrain, statistics(idx).predictTrain, ...
%             statistics(idx).scoredTest, statistics(idx).predictTest, CM_SAVE_DIR)
%     end
% end

set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying
save_stats
save(strcat(CM_SAVE_DIR, filesep, OUTPUT_STATS_FILENAME), 'save_stats');

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