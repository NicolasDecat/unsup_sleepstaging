%% Cross-validation code
% Classification algorithm: K-means clustering + Nearest centroid
% classifier (similat to k-NN classifer)
% #####################################################################
% Store output of each trial(Nf loop) into a struct containing all
% information
configuration_settings;

%% Obtain epochID using epochCounter() function
whichData = WHICH_DATA;
annotation = load(ANSWER_FILE);
label = annotation.sleepstage;
stgID = epochCounter(whichData,label);
% Training
trainingProportion = TRAINING_PERCENTAGE;
nIterations = 20;
%% Multiple iteration of randomisation and cross-validation
% Initialise result struct
block(nIterations) = struct();

for Nf = 1:nIterations  
    % Sample epoch ID from each stage and divide into training and test
    [block(Nf).trainTS,block(Nf).testTS]=epochSelect(stgID,trainingProportion);
    % trainTS and testTS are the time segment ID for training set and test set respectively.

    % The following turn the nxm matrix to 1x(n*m) matrix
    trainTS = block(Nf).trainTS.';
    trainTS = trainTS(:).';
    testTS = block(Nf).testTS.';
    testTS = testTS(:).';
    
    %% Select data of wanted time ID
    % Features used for clustering are specified in selectdata.m (passing
    % on hctsa_ops variable)
    % Timeseries used for crossval are specified here.
    for c=1:NUM_CHANNELS_USED_FOR_CROSSVAL
        if c==1
            trainMat = hctsa_ops(trainTS',:);
            testMat = hctsa_ops(testTS',:);
        else
            % Other channels (assumes the index is increment of c)
            trainMat = [trainMat hctsa_ops((trainTS*c)',:)];
            testMat = [testMat hctsa_ops((testTS*c)',:)];
        end
    end
    
    %% Clustering using training dataset
    % Record cluster ID and centre of each cluster
    n_clust = 5;
    [clustID,block(Nf).Kcentre,~] = kmeans(trainMat,n_clust,'Distance','sqeuclidean',...
                        'Display','final','Replicates',50,'MaxIter',500);
    % K-means clustering                  

    %% Classification of test dataset(Nearest centroid classifier)
    % Minimum Euclidean distance from centre/mean features of the cluster
    
    for n=1:length(testTS)
        % Calculate Euclidean distance from datapoint to each centre
        distance = sqrt(sum((bsxfun(@minus,block(Nf).Kcentre,testMat(n,:))).^2,2));
        % Find the cluster of minimum distance
        [~,clustTest(n)] = min(distance);
    end
    clustTest = clustTest';
    % Test set are classified into 5 classes(clusters) obtained from
    % clustering of training set.
    %% Find equivalent cluster-sleep stage pair
    % Method 1: Highest count = Equivalent cluster
    % Count number of datapoints from each class that are assigned to the clusters
    stgLabel = [0,1,2,3,5];
    % Initialise cluster-stage counter
    pro_no = zeros(length(stgID.nStg),length(unique(clustID)));
    for m = 1:length(unique(clustID))
        % epochID that belongs to this cluster
        % pro_id = ismember(stgID.selectID,block(Nf).trainTS(clustID==m));
        pro_id = find(clustID == m); %?
        % Actual labelled stage
        actualStg = label(trainTS(pro_id));
        %actualStg =  stgID.selectLabel(pro_id);
        for n = 1:length(stgID.nStg)
            pro_no(m,n) = sum(actualStg==stgLabel(n));
        end
        [~,equi_class(m)]=max(pro_no(m,:));
    end
    equi_stage = stgLabel(equi_class);
    %% Case: Counts are equal -> repeating equivalent class
    % if unique(equi_class)~=[1:5]  
    % pro_no?
    % end
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
    scoredTrain(Nf,:) = label(trainTS)';
    scoredTest(Nf,:) = label(testTS)';
    predictTrain(Nf,:) = block(Nf).equi_train;
    predictTest(Nf,:) = block(Nf).equi_test;
end % End Nf-th randomisation

%% Average output
% For comparing different k (changing features used as a condition)
PTrainCorrectList = [block(:).P_trainCorrect];
PTestCorrectList = [block(:).P_testCorrect];

Output(k).trainCorrect = mean(PTrainCorrectList); 
Output(k).testCorrect = mean(PTestCorrectList);
% Output(k).testError = mean(Perror);
% Output(k).trainError = mean(trainPerror);
% clearvars -except Output datamat feat_id features k complexity

%% Confusion matrix - for both train and test
% For one run(k)... can change which data to use later
%addpath('/Users/sleeping/Documents/MATLAB/unsup_sleep_staging/HCTSA/PeripheryFunctions/BF_ToBinaryClass.m')

%% Confusion matrix of train data
% Reshape scored and predict matrix
g_labelTrain = reshape(scoredTrain,1,[]);
g_clustTrain = reshape(predictTrain,1,[]);

% Labelled - make non-zero stage
g_labelTrain = g_labelTrain+1;

% Clustered - Use final clustering output
g_clustTrain = g_clustTrain+1;

% Cluster 6 becomes 5
g_labelTrain(g_labelTrain==6) = 5;
g_clustTrain(g_clustTrain==6) = 5;

% Visualise confusion matrix
figure;
plotconfusion_custom(g_labelTrain, g_clustTrain, 'Confusion Matrix - Training');
saveas(gcf, strcat(CM_SAVE_DIR, filesep, 'CM_TRN_', int2str(k), '.png'));

%% Confusion matrix of test data
% Reshape scored and predict matrix
g_labelTest = reshape(scoredTest,1,[]);
g_clustTest = reshape(predictTest,1,[]);

% Labelled - make non-zero stage
g_labelTest = g_labelTest+1;

% Clustered - Use final clustering output
g_clustTest = g_clustTest+1;

% Cluster 6 becomes 5
g_labelTest(g_labelTest==6) = 5;
g_clustTest(g_clustTest==6) = 5;

% Visualise confusion matrix
plotconfusion_custom(g_labelTest, g_clustTest, 'Confusion Matrix - Testing');
saveas(gcf, strcat(CM_SAVE_DIR, filesep, 'CM_TST_', int2str(k), '.png'));

%% Clear variables for the next run
clearvars -except Output datamat feat_id features k complexity CM_SAVE_DIR
