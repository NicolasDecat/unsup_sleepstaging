% Changes from original script
% - Changed configuration settings
% - epochSelect (instead of epochSelectFunction)
% - for ch=1:NUM_CHANNELS_TO_RUN, instead of ch=1:number_of_channels_used
% - clear max (because was used as a variable beforehand)

%% Function: Count the number of epochs in each stages and recore the epochIDs

function [statsOut testMat scoredTest predictTest Nf testTS Mean_AUC AUC_per_feature] = cross_validation_selectivefeatures(experiment, hctsa_ops, cm_save_dir, number_of_channels_used, epochSelectFunction, selective_feature,sub,v,col,FF)

%% Cross-validation code

configuration_settings;


%% Obtain epochID using epochCounter() function
whichData = WHICH_DATA;
annotation = load(ANSWER_FILE);
label = annotation.sleepstage;
stgID = epochCounter(whichData,label);
stgLab = {'W','N1','N2','N3','R'};

% Training
trainingProportion = TRAINING_PERCENTAGE;
nIterations = 10;

%% Multiple iteration of randomisation and cross-validation
% Initialise result struct
block(nIterations) = struct();
stats = struct();

for Nf = 1:nIterations
    
    
    [block(Nf).trainTS,block(Nf).testTS]=epochSelect(stgID,trainingProportion);
    
    % The following turn the nxm matrix to 1x(n*m) matrix
    trainTS = block(Nf).trainTS.';
    trainTS = trainTS(:).';
    testTS = block(Nf).testTS.';
    testTS = testTS(:).';
  
    
    %% Select data of wanted time ID
   
  
    NUM_CHANNELS_TO_RUN = v;
    
    
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

    %% SUPERVISED (same data)
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
        distance = sqrt(sum((bsxfun(@minus,block(Nf).Kcentre,testMat(n,:))).^2,2)); %%%% Calculates distance between each stage data point (Test mat) 
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
    [~, row_index] = sort(max(pro_no'), 'd');
    remaining_cluster_to_allocate=unique(clustID)';
    for s = 1:length(row_index)
        row = row_index(s);
        a_row = pro_no(row, :);
        
        while sum(a_row) > -5
            max_value = max(a_row);
            max_index = find(a_row == max_value);
            max_index = max_index(randperm(length(max_index), 1));
            cluster_allocate_index = find(ismember(remaining_cluster_to_allocate, max_index));
            if cluster_allocate_index > 0
                % The cluster ID has not been allocated yet.
                equi_class(row) = max_index;
                remaining_cluster_to_allocate(cluster_allocate_index) = [];
                break;
            else
                a_row(max_index) = -1; % Remove the maximum value because the cluster already assigned.
            end
            
        end
    end
    
    equi_stage = stgLabel(equi_class);

    
    %% Convert cluster ID into equivalent stage
    % Training
    for i=1:length(clustID)
        block(Nf).equi_train(i) = equi_stage(clustID(i));
    end
    % Test
    for j=1:length(clustTest)
        block(Nf).equi_test(j) = equi_stage(clustTest(j));
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

%% Clear variables for the next run
%clearvars -except Output datamat feat_id features k complexity CM_SAVE_DIR exps statsOut

%% Compute type1AUC

   run('type1auc_eachfeature.m')


end
