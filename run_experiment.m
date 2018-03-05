% Select features from reduced_ops.txt, make another .mat file to work on
% cross-validation.
clear all; clc;

configuration_settings;
addpath(genpath('./feature_selection'));

homedir = pwd;
% Start HCTSA tools
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

annotation = load(ANSWER_FILE);
sleep_answers = annotation.sleepstage;
% 
%% Perform feature selection (experimental)

% Include dependencies
% addpath(strcat(FSLIB_TOOLBOX_DIR, filesep, 'lib')); % dependencies
% addpath(strcat(FSLIB_TOOLBOX_DIR, filesep, 'methods')); % FS methods
% addpath(genpath(strcat(FSLIB_TOOLBOX_DIR, filesep, 'lib/drtoolbox')));
% 
% features_channel1 = datamat(334:1374, :);
% features_channel2 = datamat((1374+334):1374*2, :);
% features_channel3 = datamat(((1374*2)+334):1374*3, :);

% X_train = [features_channel1; features_channel2; features_channel3];
% Y_train = [label(334:1374); label(334:1374); label(334:1374)];

% X_train = features_channel1;
% Y_train = label(334:1374);
% numF = size(X_train, 2);

%[ranking, w] = reliefF(x_data, y_data, 20);
%[ranking, w, subset] = ILFS_auto(x_data, y_data, 4, 0 )
%ranking = mRMR(x_data, y_data, size(x_data, 2));
%[ ranking , w] = mutInfFS( X_train, Y_train, size(X_train , 2));
%[ ranking , w] = fsvFS( X_train, Y_train, size(X_train , 2) );

%Laplacian
% W = dist(X_train');
% W = -W./max(max(W)); % it's a similarity
% [lscores] = LaplacianScore(X_train, W);
% [junk, ranking] = sort(-lscocd cd gires);

% MCFS: Unsupervised Feature Selection for Multi-Cluster Data
% options = [];
% options.k = 5; %For unsupervised feature selection, you should tune
% %this parameter k, the default k is 5.
% options.nUseEigenfunction = 4;  %You should tune this parameter.
% [FeaIndex,~] = MCFS_p(X_train,numF,options);
% ranking = FeaIndex{1};
        
% ranking = spider_wrapper(X_train,Y_train,numF,lower('rfe'));
% ranking = spider_wrapper(X_train,Y_train,numF,lower('10'));
% ranking = spider_wrapper(X_train,Y_train,numF,lower('fisher'));

% Infinite Feature Selection 2015 updated 2016
% alpha = 0.5;    % default, it should be cross-validated.
% sup = 1;        % Supervised or Not
% [ranking, w] = infFS( X_train , Y_train, alpha , sup , 0 );    

% This is matlab feature selection
%[ranked, weight] = relieff(features_channel1, label(trainTS), 10, 'method', 'classification', 'categoricalx', 'on');

% Features Selection via Eigenvector Centrality 2016
% alpha = 0.5; % default, it should be cross-validated.
% ranking = ECFS( X_train, Y_train, alpha )  ;

% Regularized Discriminative Feature Selection for Unsupervised Learning
% nClass = 2;
% ranking = UDFS(X_train , nClass ); 

% BASELINE - Sort features according to pairwise correlations
% ranking = cfs(X_train);     
%         
% [B, I] = sort(ranking);
% feat_id = I;

%% Run cross-validation code
% Change the number of operations
set(0,'DefaultFigureVisible','off') % Remove this to disable the figure displaying (sometimes it could be lots of figures!)
%for k = 1:10 % k is the condition to select operation
exps = EXPS_TO_RUN; % This allow us to selectively choose which experiment to run

exp_configuration = readtable('experiment_runs.csv');
actual_exp_run = exp_configuration; 
actual_exp_run(:,:) = []; % delete all rows to accumulate the actual run

%% Create a directory to store all the result
current_dir = pwd;
output_folder = strcat(pwd, filesep, sprintf('ExpRunResult_%s',datestr(clock,'yyyymmdd_HHMMSS')));
mkdir(output_folder);

myCluster = parcluster('local')
myCluster.NumWorkers = THREAD_POOL
parpool(THREAD_POOL);

try
    parfor indx = 1:length(exps)
        exp_count = exps(indx);
        exp = exp_configuration(exp_count, :);
        
        disp(strcat('Running experiment ', int2str(exp.id), ': ', exp.name, '...'));

        max_stats = [];
        for repeat = 1:exp.repeat
            fs_algorithm = char(exp.fs_algorithm);

            selected_features = [];
            selected_feature_indexes = [];

            if (strcmp(char(exp.fs_algorithm), 'BEN') == 1)
                [selected_features, selected_feature_indexes] = fs_htsca(exp, features, feat_id);
            elseif (strcmp(char(exp.fs_algorithm), 'RELIEFF') == 1)
                [selected_features, selected_feature_indexes] = fs_reliefF(exp, datamat, sleep_answers);
            end

            % Safeguard to ensure the indexes length are always <= features
            if (length(selected_feature_indexes) > length(selected_features))
                disp(strcat('WARNING: Feature length is greater than features.', 'Experiment ', num2str(exp.id)));
                selected_feature_indexes = selected_feature_indexes(1:length(selected_features));
            end

            crossval_features = selected_features(selected_feature_indexes);
            hctsa_ops = datamat(:, crossval_features);
            statsOut = cross_validation(exp.id, hctsa_ops, output_folder);
            statsOut.complexity = exp.id;
            statsOut.features = crossval_features;

            if isempty(max_stats)
               max_stats = statsOut;
            elseif statsOut.output.trainCorrect > max_stats.output.trainCorrect
               max_stats = statsOut; 
            end
        end

        % Collect only the statistics for the maximum trainCorrect for all
        % repeats.
        statistics(indx, :) = max_stats;
        
        exp.trainCorrect = max_stats.output.trainCorrect;
        exp.testCorrect = max_stats.output.testCorrect;
        feature_str = sprintf('%.0f,' , max_stats.features);
        feature_str = feature_str(1:end-1);% strip final comma
        
        exp.feature_list{1} = feature_str;
        exp_runs(indx, :) = exp;
    end
catch ME
    delete(gcp('nocreate'));
    rethrow(ME);
end

delete(gcp('nocreate'));

%% Draw the confusion matrix for the repeat that has maximum trainCorrect
if PLOT_CONFUSION_MATRIX
    for idx = 1:length(exps)
        exp = exp_configuration(exps(idx), :);
        plot_confusion_matrix(exp.id, statistics(idx).scoredTrain, statistics(idx).predictTrain, ...
            statistics(idx).scoredTest, statistics(idx).predictTest, output_folder)
    end
end

%% Start the base line
% hctsa_ops = datamat(:,feat_id);
% statsOut = cross_validation('Baseline', hctsa_ops, LABEL_ASSOC_ALGORITHM);
% 
% disp(strcat('Baseline', ': Training percentage:',' ', num2str(statsOut.output.trainCorrect) , ...
%        ' Testing percentage:', ' ', num2str(statsOut.output.testCorrect)));
% 
% parpool(10);
% indexes = FILTERED_FEATURE_IDS; %Relief F
% thshd = EVAL_THRESHOLD;
% for top_n = NUMBER_OF_FEATURES_CHOSEN
%     disp(strcat('Running top-',num2str(top_n),'...'));
%     parfor repeats = 1:10
%         hctsa_ops = datamat(:,indexes(1:top_n));
%         stats(repeats,:) = cross_validation(strcat('Top', top_n), hctsa_ops, LABEL_ASSOC_ALGORITHM);
%     end
%     
%     parfor repeats = 1:10
%         output = stats(repeats,:).output;
%         if (output.trainCorrect - statsOut.output.trainCorrect > thshd && ...
%         output.testCorrect - statsOut.output.testCorrect > thshd)
%             disp(strcat('Top-', num2str(top_n), ' is better than baseline', ': Training percentage:',' ', num2str(output.trainCorrect) , ...
%                 ' Testing percentage:', ' ', num2str(output.testCorrect)));
%             allOneString = sprintf('%.0f,' , indexes(1:top_n));
%             allOneString = allOneString(1:end-1);% strip final comma
%             disp(allOneString);
%         end
%     end
% end
% 
% for rnd_n = NUMBER_OF_FEATURES_CHOSEN
%     disp(strcat('Running random ',num2str(rnd_n),'...'));
%     for i = 1:10
%         parfor repeats = 1:10
%             rand_idx(repeats,:) = randperm(length(indexes),rnd_n);
%             hctsa_ops = datamat(:,indexes(rand_idx(repeats,:)));
%             stats(repeats,:) = cross_validation(strcat('Random', rnd_n), hctsa_ops, LABEL_ASSOC_ALGORITHM);
%         end
% 
%         parfor repeats = 1:10
%             output = stats(repeats,:).output;
%             if (output.trainCorrect - statsOut.output.trainCorrect > thshd && ...
%                 output.testCorrect - statsOut.output.testCorrect > thshd)
%                 disp(strcat('Random ', num2str(rnd_n), ' is better than baseline', ': Training percentage:',' ', num2str(output.trainCorrect) , ...
%                     ' Testing percentage:', ' ', num2str(output.testCorrect)));
%                 allOneString = sprintf('%.0f,' , indexes(rand_idx(repeats,:)));
%                 allOneString = allOneString(1:end-1);% strip final comma
%                 disp(allOneString);
%             end
%         end
%         clear rand_idx stats;
%     end
% end
% 
set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying


    %% DBSCAN
%     allTrainMat = allTrainMat(1:1000,:);
%     
%     addpath('/Users/Zhao/Documents/MATLAB/Add-Ons/Collections/DBSCAN Clustering Algorithm/code/YPML110 DBSCAN Clustering/DBSCAN Clustering'); % dependencies
%     max_cluster = 0;
%     max_idx = [];
%     max_i = 0;
%     max_e = 0;
%     for i=2:198
%         X=allTrainMat(:,1:i);
%         epsilon = 0.1;
%         
%         while epsilon <= 1.0
%             MinPts=i+1;
%             IDX=DBSCAN(X,epsilon,MinPts);
% 
%             unique_cluster = length(unique(IDX));
%             if unique_cluster > max_cluster
%                 max_cluster = unique_cluster;
%                 max_idx = IDX;
%                 max_i = i;
%                 max_e = epsilon;
%             end
%             
%             epsilon = epsilon + 0.1;
%         end
%     end

%     disp("At " + max_i + ", max cluster " + max_cluster + " discovered. e = " + max_e);
%% Plot Results
% PlotClusterinResult(X, IDX);
% title(['DBSCAN Clustering (\epsilon = ' num2str(epsilon) ', MinPts = ' num2str(MinPts) ')']);    
% shg;


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
    plot(complexity,accuracy_train,complexity,accuracy_test)
    legend('Training','Test')
    ylabel('Accuracy [0-1]')
    xlabel('Experiment Id')
    set(gca,'XTick',(min(exps):1:max(exps)))

    saveas(gcf, strcat(output_folder, filesep, 'ACCURACY_REPORT.png'));
end

for i = 1:height(exp_runs)
    if isempty(actual_exp_run)
        actual_exp_run = exp_runs(i,:);
    else
        actual_exp_run = [actual_exp_run; exp_runs(i,:)];
    end
end

%% Sort the actual run table based on the trainCorrect descending
col_index = find(strcmpi(actual_exp_run.Properties.VariableNames,'trainCorrect'));
[sorted_table, indexes] = sortrows(actual_exp_run, col_index, 'descend');

% Print out the top 3
top_n = 3;
if top_n > height(sorted_table)
   top_n = height(sorted_table);
end

disp(sorted_table(1:top_n, :));
writetable(sorted_table, strcat(output_folder, filesep, 'experiment_summary.csv'));

