%% Function: Count the number of epochs in each stages and recore the epochIDs
% Input: - whichData // which dataset to be used ([1,5,7,13,14] for now)
%        - LabeledStage // labelled sleep stage from annotation
% Output : stgID // randomised order of epoch IDs

function statsOut = cross_validation(experiment, hctsa_ops, cm_save_dir)

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
stats = struct();

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
                        'Display','off','Replicates',50,'MaxIter',500);
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
    
    %% Match up the nearest neighbours with cluster ID
%     matched_stages = [];
%     for m = 1:length(unique(clustID))
%        clust_point = block(Nf).Kcentre(m,:);
%        
%        % Find the top 5 closest neighbours to the centroid.
%        idx = knnsearch(trainMat, clust_point, 'k', 1)';
%        
%        % Original labels
%        orig_labels = label(trainTS(idx));
       
       % Find out if any of the labels has already been matched and remove
       % them
%        indexes = find(ismember(orig_labels, matched_stages));
%        orig_labels(indexes) = [];
%        
%        % To avoid NaN (i.e. all orig_labels have been removed)
%        if (length(orig_labels) == 0)
%           orig_labels = label(trainTS(idx)); 
%        end
       
%        stage_label_most_occurence = mode(orig_labels);
%        matched_stages = [matched_stages, stage_label_most_occurence];
%        % Determine the most occurrence of the labels and associate with the
%        % cluster
%        %disp("Cluster " + m + " associated with stage " + stage_label_most_occurence);
%     end
    
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

%     addpath('/Users/Zhao/Documents/MATLAB/Add-Ons/Functions/fm_index( X, Y )/code');
%     addpath('/Users/Zhao/Documents/MATLAB/Add-Ons/Functions/The Adjusted Mutual Information/code');
    
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
end % End Nf-th randomisation

fm_train = [];
fm_test = [];
ami_train = [];
ami_test = [];
% for i=1:size(scoredTrain, 1)
%     fm_train(i) = fm_index(scoredTrain(i,:)+1, predictTrain(i, :)+1);
%     fm_test(i) = fm_index(scoredTest(i,:)+1, predictTest(i, :)+1);
%     
%     % AMI requires non-zero values
%     ami_train(i) = ami(scoredTrain(i,:)+1, predictTrain(i, :)+1);
%     ami_test(i) = ami(scoredTest(i,:)+1, predictTest(i,:)+1);
% end
% 
% stats.fm_train_score = mean(fm_train);
% stats.fm_test_score = mean(fm_test);
% stats.ami_train_score = mean(ami_train);
% stats.ami_test_score = mean(ami_test);

% disp("Mean train Fowlkes-Mallows scores for k = " + experiment + " is " + mean(fm_train));
% disp("Mean test Fowlkes-Mallows scores for k = " + experiment + " is " + mean(fm_test));
% disp("Mean train Adjusted Mutual Information scores for k = " + experiment + " is " + mean(ami_train));
% disp("Mean test Adjusted Mutual Information scores for k = " + experiment + " is " + mean(ami_test));

%% Average output
% For comparing different k (changing features used as a condition)
PTrainCorrectList = [block(:).P_trainCorrect];
PTestCorrectList = [block(:).P_testCorrect];

stats.output.trainCorrect = mean(PTrainCorrectList); 
stats.output.testCorrect = mean(PTestCorrectList);

statsOut = stats;

%% Clear variables for the next run
%clearvars -except Output datamat feat_id features k complexity CM_SAVE_DIR exps statsOut

end
