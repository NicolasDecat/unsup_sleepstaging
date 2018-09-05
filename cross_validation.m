%% Function: Count the number of epochs in each stages and recore the epochIDs
% Input: - whichData // which dataset to be used ([1,5,7,13,14] for now)
%        - LabeledStage // labelled sleep stage from annotation
% Output : stgID // randomised order of epoch IDs

function statsOut = cross_validation(experiment, hctsa_ops, cm_save_dir, number_of_channels_used, total_channels, epochSelectFunction)

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
stgLab = {'W','N1','N2','N3','R'};

s=find(label~=0); epochEnd = s(end)+2; if (epochEnd > length(label)) epochEnd = length(label); end; strcat(num2str(s(1)-1),'-',num2str(epochEnd))

% Safeguard against mismatch of annotation data and the actual epoch
% numbers.
assert(size(label, 1) == size(hctsa_ops, 1)/total_channels, "The length of the annotation must be the same as one channel length");


%% Evaluate the clustering as a whole
clust = zeros(size(hctsa_ops,1), 10);
for i = 1:10
    clust(:, i) = kmeans(hctsa_ops, i,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);
end

s=evalclusters(hctsa_ops, clust, 'silhouette');
c=evalclusters(hctsa_ops, clust, 'CalinskiHarabasz');
db=evalclusters(hctsa_ops, clust, 'DaviesBouldin');

fprintf("Full dataset: Silhouette: %.04f (%d) CalinskiHarabasz: %.04f (%d) DaviesBouldin: %.04f (%d)\n\n", ...
    s.CriterionValues(s.OptimalK), ...
    s.OptimalK, ...
    c.CriterionValues(c.OptimalK), ...
    c.OptimalK, ...
    db.CriterionValues(db.OptimalK), ...
    db.OptimalK);
        

%%

% Training
trainingProportion = TRAINING_PERCENTAGE;
nIterations = CROSSVAL_ITERATION;

%% Multiple iteration of randomisation and cross-validation
% Initialise result struct
block(nIterations) = struct();
stats = struct();

for Nf = 1:nIterations
    % Sample epoch ID from each stage and divide into training and test
    [block(Nf).trainTS,block(Nf).testTS]=epochSelectFunction(stgID,trainingProportion);
    % trainTS and testTS are the time segment ID for training set and test set respectively.

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
        
%         for stage=1:size(block(Nf).trainTS,1)
%             stage_folder = strcat(debug_folder, filesep, stgLab{stage});
%             mkdir(stage_folder);
%             
%             stage_data = block(Nf).trainTS(stage, :);
%             for k=1:length(stage_data)
%                 imageFile = strcat('ccshs_', sprintf('%03d', whichData), '_', sprintf('%04d', stage_data(k)), '.png');
%                 copyfile(strcat(DEBUG_CROSSVALIDATION_IMAGEDIR, filesep, imageFile), strcat(stage_folder, filesep, imageFile));
%             end
%         end

    end
    

    %% Select data of wanted time ID
    % Features used for clustering are specified in selectdata.m (passing
    % on hctsa_ops variable)
    % Timeseries used for crossval are specified here.
    for ch=1:number_of_channels_used
        if ch==1
            trainMat = hctsa_ops(trainTS',:);
            testMat = hctsa_ops(testTS',:);
        else
            increment=(size(hctsa_ops,1)/total_channels)*(ch-1);
            trainMat = [trainMat hctsa_ops((trainTS+increment)',:)];
            testMat = [testMat hctsa_ops((testTS+increment)',:)];
        end
    end
    
    %% Clustering using training dataset
    % Record cluster ID and centre of each cluster

    %% UNSUPERVISED
    [clustID,block(Nf).Kcentre,~] = kmeans(trainMat,NUM_CLUSTERS,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);
                   
    %% Find the max cluster that has highest silhoutte values
    clust = zeros(size(trainMat,1), 10);
    for i = 1:10
        clust(:, i) = kmeans(trainMat, i,'Distance','sqeuclidean',...
                            'Display','off','Replicates',50,'MaxIter',500);
    end
                    
    s=evalclusters(trainMat, kmeans_func, 'silhouette');
    c=evalclusters(trainMat, kmeans_func, 'CalinskiHarabasz');
    db=evalclusters(trainMat, kmeans_func, 'DaviesBouldin');
    
    total_epochs = size(trainMat, 1) + size(testMat, 1);
    
    fprintf("(Balanced N=%d) Silhouette: %.04f (%d) CalinskiHarabasz: %.04f (%d) DaviesBouldin: %.04f (%d)\n\n", ...
        total_epochs, ...
        s.CriterionValues(s.OptimalK), ...
        s.OptimalK, ...
        c.CriterionValues(c.OptimalK), ...
        c.OptimalK, ...
        db.CriterionValues(db.OptimalK), ...
        db.OptimalK);
    
    %% SUPERVISED (same data)
    template = templateSVM(...
    'KernelFunction', 'linear', ...
    'PolynomialOrder', [], ...
    'KernelScale', 'auto', ...
    'BoxConstraint', 1, ...
    'Standardize', true);

     SVMModel = fitcecoc(trainMat, label(trainTS),'Learners', template, ...
         'Coding', 'onevsone', 'ClassNames', [0; 1; 2; 3; 5]);
     svmTrain = predict(SVMModel, trainMat);
     svmTest = predict(SVMModel, testMat);
     
    %% Classification of test dataset(Nearest centroid classifier)
    % Minimum Euclidean distance from centre/mean features of the cluster
    
    for n=1:length(testTS)
        % Calculate Euclidean distance from datapoint to each centre
        distance = sqrt(sum((bsxfun(@minus,block(Nf).Kcentre,testMat(n,:))).^2,2));
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
        % epochID that belongs to this cluster
        % pro_id = ismember(stgID.selectID,block(Nf).trainTS(clustID==m));
        pro_id = find(clustID == m); %?
        % Actual labelled stage
        actualStg = label(trainTS(pro_id));
        %actualStg =  stgID.selectLabel(pro_id);
        for n = 1:length(stgID.nStg)
            pro_no(m,n) = sum(actualStg==stgLabel(n));
        end
        %[~,equi_class(m)]=max(pro_no(m,:));
    end
    
    % Now we give the stage that has the highest sum choose first.
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
        
        %% Find the top X closest neighbours to each state.
%         for m = 1:length(unique(clustID))
%             stage_name = print_stage_name{m};
%             clust_point = block(Nf).Kcentre(m,:);
%             idx = knnsearch(trainMat, clust_point, 'k', 5)';
%             ts_idx = trainTS(idx);
%             
%             stage_folder = strcat(debug_folder, filesep, stage_name);
%             
%             ll=label(ts_idx);
%             ll = ll+1;
%             ll(ll == 6) = 5;
%             
%             for s = 1:length(ts_idx)
%                 imageFile = strcat('ccshs_', sprintf('%03d', whichData), '_', sprintf('%04d', ts_idx(s)), '.png');
%                 
%                 suffix = '_CORRECT';
%                 if (stgLab{ll(s)} ~= stage_name)
%                     suffix = strcat('_INCORRECT_', stgLab{ll(s)});
%                 end
%                 
%                 imageToFile = strcat('ccshs_', sprintf('%03d', whichData), '_', sprintf('%04d', ts_idx(s)), '_CLOSEST', suffix, '.png');
%                 copyfile(strcat(DEBUG_CROSSVALIDATION_IMAGEDIR, filesep, imageFile), strcat(stage_folder, filesep, imageToFile));
%             end
%         end
    end
    
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

    %fprintf('trainCorrect: %0.2f testCorrect: %0.2f\n', block(Nf).P_trainCorrect, block(Nf).P_testCorrect);
    
    %% Confusion matrix input
    stats.scoredTrain(Nf,:) = label(trainTS)';
    stats.scoredTest(Nf,:) = label(testTS)';
    stats.predictTrain(Nf,:) = block(Nf).equi_train;
    stats.predictTest(Nf,:) = block(Nf).equi_test;
    stats.svmPredictTrain(Nf, :) = svmTrain';
    stats.svmPredictTest(Nf, :) = svmTest';
    
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

%% Clear variables for the next run
%clearvars -except Output datamat feat_id features k complexity CM_SAVE_DIR exps statsOut

end
