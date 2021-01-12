%%%%%
%%%%%

%% classification performance for each feature: kmeans for each binary classifier


% Find how many epochs for each sleep stage
stgL = size(testMat{1},1);
stgL = stgL/5;          

% Get epochs indices for each stage
wake = (1:stgL);
N1 = ((stgL+1):(stgL*2));
N2 = ((stgL*2+1):(stgL*3));
N3 = ((stgL*3+1):(stgL*4));
rem = ((stgL*4+1):(stgL*5));

% All binary classifiers
allpairs = [{[wake;N1]} {[wake;N2]} {[wake;N3]} {[wake;rem]} {[N1;N2]} {[N1;N3]} {[N1;rem]} {[N2;N3]} {[N2;rem]} {[N3;rem]}];

for Nit = 1:length(testMat)

    for C = 1:10    % For each binary classifier

            % Choose one classifier
            classifier = allpairs{C};     

            for test = 1:stgL           % each epoch will be taken out and used as testing epoch

            % Get epoch indices for both stages of the classifier
            stage1 = classifier(1,:);    
            stage2 = classifier(2,:); 

            stages = [ones(1,stgL) 2*ones(1,stgL)];  % vector of 1 and 2 (used later for kmeans)

            % Testing data (1 epoch for each stage is taken out: "leave-1-out" strategy for cross validation)
            test_stage1 = stage1(test);   % The only piece of information that we change for each 'test' iteration
            test_stage2 = stage2(test);

            % Training data (N-1 epochs for each stage)
            stage1 = setdiff(stage1,test_stage1);   % Leave one epoch out ("test" epoch)
            stage2 = setdiff(stage2,test_stage2);  


            %% kmeans
            NUMCLUST = 2;

            % hctsa response of training
            hctsa_resp1 = testMat{Nit}(stage1,:);   % hctsa resp for all 10 epochs of stage 1, including all features
            hctsa_resp2 = testMat{Nit}(stage2,:);

            hctsa_resp = [hctsa_resp1;hctsa_resp2];

            % hctsa response of testing
            hctsa_test1 = testMat{Nit}(test_stage1,:);
            hctsa_test2 = testMat{Nit}(test_stage2,:);

            hctsa_test = [hctsa_test1;hctsa_test2];

            [clustID,blokk(test).Kcentre] = kmeans(hctsa_resp,NUMCLUST,'Distance','sqeuclidean',...
                        'Display','off','Replicates',50,'MaxIter',500);


            %% Classification of test dataset (Nearest centroid classifier)

            for n=1:size(hctsa_test,1)   % for both the testing epoch of stage 1 and of stage 2
                distance = sqrt(sum((bsxfun(@minus,blokk(test).Kcentre,hctsa_test(n,:))).^2,2));
                [m,clusttest(n)] = min(distance);   % clusttest = cluster ID assigned to testing epoch of stage 1 and 2 based on distance between Kcenters and testing epoch hctsa response
            end

            CLUSTTEST(test,1:2) = clusttest;

            %% Find equivalent cluster-sleep stage pair

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


            % To know which of Clust ID 1 and 2 corresponds to stage 1 and 2
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
          
    end
    
end


disp('ok')

  


