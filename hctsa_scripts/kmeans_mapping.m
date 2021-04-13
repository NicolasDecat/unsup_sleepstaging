% Changes from original script
% - Changed configuration settings
% - epochSelect (instead of epochSelectFunction)
% - for ch=1:NUM_CHANNELS_TO_RUN, instead of ch=1:number_of_channels_used
% - clear max (because was used as a variable beforehand)

%% Function: Count the number of epochs in each stages and recore the epochIDs

function [statsOut testMat scoredTest predictTest Nf Iteration NumChannels Dataset Sleep_stage Testing_accuracy AUC testTS] = cross_validation_selectivefeatures(experiment, hctsa_ops, cm_save_dir, number_of_channels_used, epochSelectFunction, selective_feature,sub,v,col)

%% Cross-validation code

configuration_settings;

% AUC parameters
origlabels = [];
clusterdecision = [];
TestingAcc = [];
AUC = [];
NumIteration = [];
stgAUC = [];
Dataset = [];
NumChannels = [];


%% Obtain epochID using epochCounter() function
whichData = WHICH_DATA;
annotation = load(ANSWER_FILE);
label = annotation.sleepstage;
stgID = epochCounter(whichData,label);
stgLab = {'W','N1','N2','N3','R'};

% Training
trainingProportion = TRAINING_PERCENTAGE;
nIterations = 1;

%% Multiple iteration of randomisation and cross-validation
% Initialise result struct
block(nIterations) = struct();
stats = struct();

for Nf = 1:nIterations
    
    disp(sprintf('Dataset %s,Channel %s, Iteration %s',sub, num2str(v),num2str(Nf)))
    
    [block(Nf).trainTS,block(Nf).testTS]=epochSelect(stgID,trainingProportion);
    
    % The following turn the nxm matrix to 1x(n*m) matrix
    trainTS = block(Nf).trainTS.';
    trainTS = trainTS(:).';
    testTS = block(Nf).testTS.';
    testTS = testTS(:).';

    if DEBUG_CROSSVALIDATION
        debug_folder =  strcat('Exp_', sprintf('%d', experiment), '_iteration_', num2str(Nf), '_channel_', num2str(number_of_channels_used));
        if exist(debug_folder)==7
           rmdir(debug_folder, 's');
        end
        mkdir(debug_folder);       
        
        for stage=1:size(block(Nf).trainTS,1)
            stage_folder = strcat(debug_folder, filesep, stgLab{stage});
            mkdir(stage_folder);
            
            stage_data = block(Nf).trainTS(stage, :);
            for k=1:length(stage_data)
                imageFile = strcat('ccshs_', sprintf('%03d', whichData), '_', sprintf('%04d', stage_data(k)), '.png');
                copyfile(strcat(DEBUG_CROSSVALIDATION_IMAGEDIR, filesep, imageFile), strcat(stage_folder, filesep, imageFile));
            end
        end

    end
  
    
    %% Select data of wanted time ID
   
    if v==1
            NUM_CHANNELS_TO_RUN = [1];
        elseif v == 2
            NUM_CHANNELS_TO_RUN = [2];
        elseif v== 3
            NUM_CHANNELS_TO_RUN = [3];
    end
    
    
    for ch=1:NUM_CHANNELS_TO_RUN

        if ch==1
            trainMat = hctsa_ops(trainTS',:);
            testMat = hctsa_ops(testTS',:);
        else
            % increment=(size(hctsa_ops,1)/3)*(ch-1);
            increment=round((size(hctsa_ops,1)/7))*(ch-1);
            trainMat = [trainMat hctsa_ops((trainTS+increment)',:)];
            testMat = [testMat hctsa_ops((testTS+increment)',:)];
        end
    end
    
    
    %% Clustering using training dataset
    % Record cluster ID and centre of each cluster

    %% UNSUPERVISED
    [clustID,block(Nf).Kcentre,sse] = kmeans(trainMat,NUM_CLUSTERS,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);


%     %% SUPERVISED (same data)
%     template = templateSVM(...
%     'KernelFunction', 'linear', ...
%     'PolynomialOrder', [], ...
%     'KernelScale', 'auto', ...
%     'BoxConstraint', 1, ...
%     'Standardize', true);
% 
%      SVMModel = fitcecoc(trainMat, label(trainTS),'Learners', template, ...
%          'Coding', 'onevsone', 'ClassNames', [0; 1; 2; 3; 5]);
%      svmTrain = predict(SVMModel, trainMat);
%      svmTest = predict(SVMModel, testMat);
     
     
    %% Classification of test dataset(Nearest centroid classifier)
    % Minimum Euclidean distance from centre/mean features of the cluster
    
    for n=1:length(testTS)
        % Calculate Euclidean distance from datapoint to each centre
        distance = sqrt(sum((bsxfun(@minus,block(Nf).Kcentre,testMat(n,:))).^2,2)); %%%% Calculate distance between each stage data point (Test mat) 
                                                                                    %%%% and center of cluster from kmeans. Data point with least distance from a center goes to cluster of that center. ClustTest = predictTest, = the cluster decisions
        % Find the cluster of minimum distance              
        [~,clustTest(n)] = min(distance);
    end
      
%% Test set are classified into 5 classes(clusters) obtained from
    % clustering of training set.
    %% Find equivalent cluster-sleep stage pair
    % Method 1: Highest count = Equivalent cluster
    % Count number of datapoints from each class that are assigned to the clusters
    stgLabel = [0,1,2,3,5];
    % Initialise cluster-stage counter
    pro_no = zeros(length(stgID.nStg),length(unique(clustID)));
    for m = 1:length(unique(clustID))
       
        pro_id = find(clustID == m); 
        % Actual labelled stage
        actualStg = label(trainTS(pro_id));
        %actualStg =  stgID.selectLabel(pro_id);
        for n = 1:length(stgID.nStg)
            pro_no(m,n) = sum(actualStg==stgLabel(n));
        end
    end
    
    % Now we give the stage that has the highest sum choose first.
    clear max
    [~, row_index] = sort(max(pro_no'), 'd');    % From clust ID that has the highest number of epochs of the most frequent stage, to least epoch
    remaining_cluster_to_allocate=unique(clustID)';   % all 5 cluster IDs
    for s = 1:length(row_index)    % For each of the 5 cluster IDs
        row = row_index(s);         % --> Start with the cluster ID that has the highest number of epoch of the most frequent stage 
        a_row = pro_no(row, :);   % Show their number of epochs for each stage
        
        while sum(a_row) > -5
            max_value = max(a_row);    % Number of epoch of the most frequent stage
            max_index = find(a_row == max_value);  % Most frequent stage (1 2 3 4 5), 4 being N3
            max_index = max_index(randperm(length(max_index), 1));  % If 2 have same high number, take rand
            cluster_allocate_index = find(ismember(remaining_cluster_to_allocate, max_index));  % cluster to allocate is a cluster that has not been allocated yet 
            if cluster_allocate_index > 0
                % The cluster ID has not been allocated yet.
                equi_class(row) = max_index;   % row =  each cluster ID is given the name of the most frequent stage (1 2 3 4 5)
                remaining_cluster_to_allocate(cluster_allocate_index) = [];
                break;
            else
                a_row(max_index) = -1; % Remove the maximum value because the cluster already assigned.
            end
            
        end
    end
    
    equi_stage = stgLabel(equi_class);   % convert cluster ID into equivalent class. 

    if DEBUG_CROSSVALIDATION
        print_stage = equi_stage+1;
        print_stage(print_stage == 6) = 5;

        print_stage_name = stgLab(print_stage);

        fileID = fopen(strcat(debug_folder, filesep, 'results.txt'),'w');
        fmt = [repmat('%4d ', 1, size(pro_no,2)-1), '%4d\n'];
        fprintf(fileID, fmt, pro_no.');   %transpose is important!
        fprintf(fileID, '\n');
        
        for l=1:length(print_stage_name)
           fprintf(fileID, 'Cluster %d => %s\n ', l, print_stage_name{l});
        end
    end
    
    
    %% Convert cluster ID into equivalent stage
    % Training
    for i=1:length(clustID)
        block(Nf).equi_train(i) = equi_stage(clustID(i));
    end
    % Test
    for j=1:length(clustTest)
        block(Nf).equi_test(j) = equi_stage(clustTest(j));
    end

    if DEBUG_CROSSVALIDATION
        Index = trainTS';
        Answers = label(trainTS);
        Predicts = block(Nf).equi_train';

        csvTable = table(Index,Answers,Predicts);
        writetable(csvTable, strcat(debug_folder, filesep, 'summary.csv'));
        
        plot_confusion_matrix(experiment, label(trainTS)', block(Nf).equi_train, ...
            label(testTS)', block(Nf).equi_test, debug_folder)
    end
    
    %% Percentage correct
    trainCorrect = sum(block(Nf).equi_train==label(trainTS)');
    testCorrect = sum(block(Nf).equi_test==label(testTS)');
    block(Nf).P_trainCorrect = trainCorrect/(length(trainTS));
    block(Nf).P_testCorrect = testCorrect/(length(testTS));
    
    %% Confusion matrix input
    stats.scoredTrain(Nf,:) = label(trainTS)';
    stats.scoredTest(Nf,:) = label(testTS)';
    stats.predictTrain(Nf,:) = block(Nf).equi_train;
    stats.predictTest(Nf,:) = block(Nf).equi_test;
%     stats.svmPredictTrain(Nf, :) = svmTrain';
%     stats.svmPredictTest(Nf, :) = svmTest';
    
    assert(size(trainMat, 2) == size(testMat, 2));
    stats.totalFeatures(Nf, :) = size(trainMat, 2);
end % End Nf-th randomisation

%% Average output
% For comparing different k (changing features used as a condition)
PTrainCorrectList = [block(:).P_trainCorrect];
PTestCorrectList = [block(:).P_testCorrect];

stats.output.trainCorrect = mean(PTrainCorrectList); 
stats.output.testCorrect = mean(PTestCorrectList);

statsOut = stats;

scoredTest = stats.scoredTest;
predictTest = stats.predictTest;
% scoredTest = stats.scoredTest;   % If want to switch to supervised learning
% predictTest = stats.svmPredictTest;


%% Clear variables for the next run
%clearvars -except Output datamat feat_id features k complexity CM_SAVE_DIR exps statsOut

%% Run AUC
   run('type1aucc.m')


end
