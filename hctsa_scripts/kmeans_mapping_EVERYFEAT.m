% Changes from original script
% - Changed configuration settings
% - epochSelect (instead of epochSelectFunction)
% - for ch=1:NUM_CHANNELS_TO_RUN, instead of ch=1:number_of_channels_used
% - clear max (because was used as a variable beforehand)

%% Function: Count the number of epochs in each stages and recore the epochIDs

function [statsOut testMat scoredTest predictTest testTS] = cross_validation_selectivefeatures(experiment, hctsa_ops, cm_save_dir, number_of_channels_used, epochSelectFunction, selective_feature,sub,v,col)

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
nIterations = 300;

%% Multiple iteration of randomisation and cross-validation
% Initialise result struct
block(nIterations) = struct();
stats = struct();

scoredTest = NaN(1,numel(annotation.sleepstage));
predictTest = NaN(1,numel(annotation.sleepstage));

for Nf = 1:nIterations
    
    disp(sprintf('Dataset %s,Channel %s, Iteration %s',sub, num2str(v),num2str(Nf)))
    
    [block(Nf).trainTS,block(Nf).testTS]=epochSelect(stgID,trainingProportion);
    
    % The following turn the nxm matrix to 1x(n*m) matrix
    trainTS = block(Nf).trainTS.';
    trainTS = trainTS(:).';
    testTS = block(Nf).testTS.';
    testTS = testTS(:).';
    
    whichTS(Nf,:) = testTS;

  
    
    %% Select data of wanted time ID
  
    
    for ch=1:v

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



     
     
    %% Classification of test dataset(Nearest centroid classifier)
    
    for n=1:length(testTS)
        distance = sqrt(sum((bsxfun(@minus,block(Nf).Kcentre,testMat(n,:))).^2,2));                                                                    
        [~,clustTest(n)] = min(distance);
    end
      
    %% Find equivalent cluster-sleep stage pair


    stgLabel = [0,1,2,3,5];
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

for i = 1:nIterations
    scoredTest(whichTS(i,:)) = stats.scoredTest(i,:);  % AASM label
    predictTest(whichTS(i,:)) = stats.predictTest(i,:);  % Cluster decisions
end



end
