%%%%%
%%%%%   Unsupervised + supervised (SVM) clustering, to compare them in
%%%%%   terms of AUC and testing_accuracy. 

%% classification performance for each feature: kmeans for each binary classifier

% Find how many testing epochs for each sleep stage
stgL = size(testMat{1},1);
stgL = stgL/5;          

% Get the indices of testing epochs for each sleep stage
wake = (1:stgL);
N1 = ((stgL+1):(stgL*2));
N2 = ((stgL*2+1):(stgL*3));
N3 = ((stgL*3+1):(stgL*4));
rem = ((stgL*4+1):(stgL*5));

% All binary classifiers for testing 
allpairs = [{[wake;N1]} {[wake;N2]} {[wake;N3]} {[wake;rem]} {[N1;N2]} {[N1;N3]} {[N1;rem]} {[N2;N3]} {[N2;rem]} {[N3;rem]}];


% Parameters for supervised clustering: Find how many training epochs for each sleep stage
stgL_train = size(trainMat{1},1);
idx_train = 1:stgL_train;  
stgL_train = stgL_train/5; 

% Get the indices of training epochs for each sleep stage
wakeT = (1:stgL_train);
N1T = ((stgL_train+1):(stgL_train*2));
N2T = ((stgL_train*2+1):(stgL_train*3));
N3T = ((stgL_train*3+1):(stgL_train*4));
remT = ((stgL_train*4+1):(stgL_train*5));

allpairs_training = [{[wakeT;N1T]} {[wakeT;N2T]} {[wakeT;N3T]} {[wakeT;remT]} {[N1T;N2T]} {[N1T;N3T]} {[N1T;remT]} {[N2T;N3T]} {[N2T;remT]} {[N3T;remT]}];

%%%% For each iteration, for each classifier and for each 'test' iteration

for Nit = 1:length(testMat)    % for each iteration

    for C = 1:10               % For each binary classifier      
            
            % Choose one classifier with testing epochs
            classifier = allpairs{C}; 
            
            % Choose one classifier with training eopchs (used for supervised clustering)
            stages_SVMtrain = allpairs_training{C};
    
            
            %% Unsupervised clustering

            for test = 1:stgL           % for each test iteration, a different epoch is chosen as testing epoch
                
            % Get epoch indices for both stages of the classifier
            stage1 = classifier(1,:);    
            stage2 = classifier(2,:); 

            stages = [ones(1,stgL) 2*ones(1,stgL)];  % vector of 1s and 2s (used later for kmeans)

            % Testing data (1 epoch for each stage is taken out: "leave-1-out" strategy for cross validation)
            test_stage1 = stage1(test);   % The only piece of information that we change for each 'test' iteration
            test_stage2 = stage2(test);

            % Training data (N-1 epochs for each stage)
            stage1 = setdiff(stage1,test_stage1);   % Leave one epoch out ("test" epoch)
            stage2 = setdiff(stage2,test_stage2);  
  
            NUMCLUST = 2;

            % hctsa response of training
            hctsa_resp1 = testMat{Nit}(stage1,:);   % hctsa resp for all 10 training epochs of stage 1, including all features
            hctsa_resp2 = testMat{Nit}(stage2,:);

            hctsa_resp = [hctsa_resp1;hctsa_resp2];

            % hctsa response of testing
            hctsa_test1 = testMat{Nit}(test_stage1,:);
            hctsa_test2 = testMat{Nit}(test_stage2,:);

            hctsa_test = [hctsa_test1;hctsa_test2];

            [clustID,blokk(test).Kcentre] = kmeans(hctsa_resp,NUMCLUST,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);


                    
            %%% Classification of test dataset (Nearest centroid classifier)

            for n=1:size(hctsa_test,1)   % for both the testing epoch of stage 1 and of stage 2
                distance = sqrt(sum((bsxfun(@minus,blokk(test).Kcentre,hctsa_test(n,:))).^2,2));
                [m,clusttest(n)] = min(distance);   % clusttest = cluster ID assigned to testing epoch of stage 1 and 2 based on distance between Kcenters and testing epoch hctsa response
            end

            CLUSTTEST(test,1:2) = clusttest;

            
            %%% Find equivalent cluster-sleep stage pair

            stgLabel = [1,2];
            pro_no = zeros(length(stgLabel),length(unique(clustID)));  % length(stg.nStG) = number of stages (5, always)

            for m = 1:length(unique(clustID))     % for each cluster ID returned by kmeans (arbitrary IDs, just group number 1, 2 , 3, 4, 5 (2, always)
                pro_id = find(clustID == m);      % find the indices of epochs that are assigned to Cluster number 1 

                actualStg = stages(pro_id);       % get the actual, original labels (1 or 2) of the training epochs for each clust ID.

                for n = 1:2   % For each of the 2 stages 
                    pro_no(m,n) = sum(actualStg==stgLabel(n));       % For each clust ID (1 2 3 4 5), find how many N1, N2, N3, W and REM epochs 
                end
            end

            CLASSEPOCH{test} = pro_no;  % Will be used to compute percentage accuracy

            % Element 1 = number of epochs from stage 1 attributed to clustID 1
            % Element 2 = number of epochs from stage 2 attributed to clustID 1
            % Element 3 = number of epochs from stage 1 attributed to clustID 2
            % Element 4 = number of epochs from stage 2 attributed to clustID 2


            %%%  To know which of Clust ID 1 and 2 corresponds to stage 1 and 2
            
            clear max
            [~, row_index] = sort(max(pro_no'), 'd');    % From clust ID that has the highest number of epochs of the most frequent stage, to least epoch
            remaining_cluster_to_allocate=unique(clustID)';   % all 2 cluster IDs
            for s = 1:length(row_index)    % For each of the 2 cluster IDs
                    row = row_index(s);         % --> Start with the cluster ID that has the highest number of epoch of the most frequent stage 
                    a_row = pro_no(row, :);   % Show their number of epochs for each stage

                    while sum(a_row) > -5
                        max_value = max(a_row);    % Number of epoch of the most frequent stage
                        max_index = find(a_row == max_value);  % Most frequent stage (1 or 2)
                        max_index = max_index(randperm(length(max_index), 1));  % If 2 have same high number, take rand
                        cluster_allocate_index = find(ismember(remaining_cluster_to_allocate, max_index));  % cluster to allocate is a cluster that has not been allocated yet 
                        if cluster_allocate_index > 0
                            % The cluster ID has not been allocated yet.
                            equi_class(row) = max_index;   % row =  each cluster ID is given the name of the most frequent stage (1 or 2)
                            remaining_cluster_to_allocate(cluster_allocate_index) = [];
                            break;
                        else
                            a_row(max_index) = -1; % Remove the maximum value because the cluster already assigned.
                        end

                    end
                end

                equi_stage(test,1:2) = stgLabel(equi_class);  

            end


            CLASSEPOCH = CLASSEPOCH';

            for E = 1:length(equi_stage)

                 [max_stage1,ii1] = max(CLASSEPOCH{E}(1,:));  % maximum value (dominant label) and index (which class) of stage 1 
                 [max_stage2,ii2] = max(CLASSEPOCH{E}(2,:));  % maximum value (dominant label) and index (which class) of stage 2
                 sum_stage1 = sum(CLASSEPOCH{E}(1,:));
                 sum_stage2 = sum(CLASSEPOCH{E}(2,:));

                 Perc_stage1 = (max_stage1/sum_stage1)*100;
                 Perc_stage2 = (max_stage2/sum_stage2)*100;

                 % Stage 1
                 if ii1 == 1 && equi_stage(E,1) == 1 || ii1 == 2 && equi_stage(E,2) == 1   % if most labels are 1 and class is 1
                     Percentage_correct_stage1(E) = Perc_stage1;
                 else                                   % if dominant label is in the opposite class
                     Percentage_correct_stage1(E) = 100 - Perc_stage1;
                 end

                 % Stage 2     
                 if ii2 == 1 && equi_stage(E,1) == 2 || ii2 == 2 && equi_stage(E,2) == 2   % if most labels are 2 and class is 2
                     Percentage_correct_stage2(E) = Perc_stage2;
                 else                                  
                     Percentage_correct_stage2(E) = 100 - Perc_stage2;
                 end

           end

           PERC_STAGE1 = mean(Percentage_correct_stage1);
           PERC_STAGE2 = mean(Percentage_correct_stage2);

           % Training accuracy
           PERC_STAGE(Nit,C) = mean([PERC_STAGE1 PERC_STAGE2]);    % rows = each iteration. Col = each classifier

           % Testing accuracy
           for j=1:length(CLUSTTEST)
               predict_test(j,1:2) = equi_stage(j,CLUSTTEST(j,:));
           end
           
           test_correct_stage1 = ((sum(predict_test(:,1)==1'))/stgL)*100;
           test_correct_stage2 = ((sum(predict_test(:,2)==2'))/stgL)*100;
           test_correct(Nit,C) = mean([test_correct_stage1 test_correct_stage2]);
          
           
     
         %% Supervised: SVM for each classifier (use the labels to create a prediction model that classifies epochs)

         % Get training epoch ID for both stages of the classifier
         stagesSVM1_train = stages_SVMtrain(1,:);
         stagesSVM2_train = stages_SVMtrain(2,:);
         
         stages_SVMtrain = [stagesSVM1_train stagesSVM2_train];
         
         % Get testing epoch indices for both stages of the classifier
         stageSVM1 = classifier(1,:);    
         stageSVM2 = classifier(2,:); 
         
         stagesSVM = [stageSVM1 stageSVM2];

         % Generate the right labels for training (1s and 2s)
         stages_train = [ones(1,stgL_train) 2*ones(1,stgL_train)];
         
         
         % hctsa response of training epochs (will be used to train algo)
         hctsa_SVMresp1 = trainMat{Nit}(stagesSVM1_train,:);  
         hctsa_SVMresp2 = trainMat{Nit}(stagesSVM2_train,:);

         hctsa_SVMresp = [hctsa_SVMresp1;hctsa_SVMresp2];   % 1:11: hctsa resp from stage1, then 12:22 from stage2

         % SVM for binary classification: first half of hctsa_SVMresp is hctsa responses
         % associated with Stage1, second half is hctsa responses
         % associated with Stage2. The model will then associate Stage1 and
         % Stage2 to their specific hctsa responses
         % 1 SVM model per classifier, because we take different training
         % epochs (depending on which classifier is used)
         SVMModel{C} = fitcsvm(hctsa_SVMresp,stages_train);  

         % Get the decisions generated by SVM model for the training and
         % testing data
         svmTrain(:,C) = predict(SVMModel{C}, trainMat{Nit}(stages_SVMtrain,:));  % how SVM classifies training epochs used to create the model
         svmTest(:,C) = predict(SVMModel{C}, testMat{Nit}(stagesSVM,:));       % how SVM classifies new epochs (testing epochs) based on training

         % Re-store these variables
         stats.svmPredictTrain(C, :) = svmTrain(:,C)';  
         stats.svmPredictTest(C, :) = svmTest(:,C)';
         
         % Get the original labels (How the SVM is supposed to label if it is perfect)
         stats.scoredTrain = stages_train;  
         stats.scoredTest = stages;   

         % Calculate the percentage of correctly labelled training and testing epochs
         iteration_svm_training_accuracy(C,Nit) = (((sum((stats.scoredTrain == stats.svmPredictTrain(C,:))'))/size(stats.scoredTrain, 2))')*100;
         iteration_svm_testing_accuracy(C,Nit) = (((sum((stats.scoredTest == stats.svmPredictTest(C,:))'))/size(stats.scoredTest, 2))')*100;

  
           
    end
    

end


% Get the mean SVM training and testing accuracy for each classifier
for Cl = 1:10   % For each of the 10 classifiers
    iteration_svm_training_accuracy_MEAN(Cl) = mean(iteration_svm_training_accuracy(Cl,:));  
    iteration_svm_testing_accuracy_MEAN(Cl) = mean(iteration_svm_testing_accuracy(Cl,:));
end
    

disp('ok')

  


