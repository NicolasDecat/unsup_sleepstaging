
%% Obtain the distinctive features for each cluster
%%% --> OVA for each cluster, ranksum an ranking according to p value


%%% 1) Get the 5 clusters: (kmeans) and the feature values
%%% 2) Do sumrank for each feature, for each of the 5 OVA classifier
%%% 3) Rank fatures according to their log(p) value

Subs = {'001','005','439','458' '596' '748' '749' '752' '604' '807' '821' '870'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};


for D = 1:length(Subs)
    
    sub = Subs{D};

    % File selection
    WHICH_DATA = str2num(sprintf('%s',sub)); 

    % go to current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
   
    % Load TS_DataMat and Operations, and convert to 5603
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ID5603')
   
    load('HCTSA_N','Operations');
    all_op = Operations;
    EQUI = find(ismember(Operations.ID,ID5603));
    Operations = Operations(EQUI,:);
    all_op = Operations;
    
    datam = load('HCTSA_N','TS_DataMat');
    datam = datam.TS_DataMat(:,EQUI);


    % How many channels
    v = 1;   

    % Load data matrix for one feature
    datamat = datam;    
    feat_id = 1:size(datamat,2);          
    hctsa_ops = datamat(:,feat_id);      


    %% Run k-means
    
    whichData = WHICH_DATA;
    annotation = load(sprintf('ccshs_1800%s_annot.mat',sub),'sleepstage');
    label = annotation.sleepstage;
    trainingProportion = 0.7;
    
    stgID = epochCounter(whichData,label);
    stgLab = {'W','N1','N2','N3','R'};
    
    nIterations = 100;
    
    block(nIterations) = struct();

    
    for Nf = 1:nIterations
        
        while true 
        
            [block(Nf).trainTS,block(Nf).testTS]=epochSelect(stgID,trainingProportion);

            trainTS = block(Nf).trainTS.';
            trainTS = trainTS(:).';
            testTS = block(Nf).testTS.';
            testTS = testTS(:).';


             for ch=1:v

                if ch==1
                    trainMat = hctsa_ops(trainTS',:);
                    testMat = hctsa_ops(testTS',:);
                else
                    increment=round((size(hctsa_ops,1)/7))*(ch-1);
                    trainMat = [trainMat hctsa_ops((trainTS+increment)',:)];
                    testMat = [testMat hctsa_ops((testTS+increment)',:)];
                end

            end

             %% k-means

             NUM_CLUSTERS = 5;

            [clustID,block(Nf).Kcentre,sse,Dist] = kmeans(trainMat,NUM_CLUSTERS,'Distance','sqeuclidean',...
                    'Display','off','Replicates',50,'MaxIter',500);


             %% Sequential maximum matching

             for n=1:length(testTS)

                % Calculate Euclidean distance from datapoint to each centre
                distance = sqrt(sum((bsxfun(@minus,block(Nf).Kcentre,testMat(n,:))).^2,2)); %%%% Calculate distance between each stage data point (Test mat) 
                                                                                            %%%% and center of cluster from kmeans. Data point with least distance from a center goes to cluster of that center. ClustTest = predictTest, = the cluster decisions
                % Find the cluster of minimum distance              
                [~,clustTest(n)] = min(distance);

             end


            %%% Find equivalent cluster-sleep stage pair

            stgLabel = [0,1,2,3,5];

            pro_no = zeros(length(stgID.nStg),length(unique(clustID)));

            for m = 1:length(unique(clustID))

                pro_id = find(clustID == m); 
                % Actual labelled stage
                actualStg = label(trainTS(pro_id));
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

            equi_stage = stgLabel(equi_class);   

            %% Now, calculate ranksum

            %%%% 1) Get cluster names and the feature values

             CLUST = equi_stage(clustTest);   % Cluster mapped to stage

             % Name of cluster
             CLUST_NAME_idx = equi_class(clustTest);
             NAME = {'CW','C1','C2','C3','CR'};
             CLUST_NAME = NAME(CLUST_NAME_idx);

             % Feature values
             for i = 1:numel(testTS)
                FeatVal(i,:)= testMat(i,:);
             end


            %%%% 2) Calculate ranksum

            % Get index of epochs for each cluster
            CW = find(contains(CLUST_NAME,'CW'));
            C1 = find(contains(CLUST_NAME,'C1'));
            C2 = find(contains(CLUST_NAME,'C2'));
            C3 = find(contains(CLUST_NAME,'C3'));
            CR = find(contains(CLUST_NAME,'CR'));
          
          % If one CLUST_NAME is empty, then redo the procedure (new
          % selection of epochs). Else, continue (break the while loop)
            if ~isempty(CW) == 1 && ~isempty(C1) == 1 && ~isempty(C2) == 1 && ~isempty(C3) == 1 && ~isempty(CR) == 1
                break
            end
        
        end
        
        % Set the 5 OVA classifiers (IDs for feature values)
        CW_rest = [{CW};{[C1 C2 C3 CR]}];
        C1_rest = [{C1};{[CW C2 C3 CR]}];
        C2_rest = [{C2};{[C1 CW C3 CR]}];
        C3_rest = [{C3};{[C1 C2 CW CR]}];
        CR_rest = [{CR};{[C1 C2 C3 CW]}];
        
        Classif = [CW_rest C1_rest C2_rest C3_rest CR_rest];
        
        for C = 1:5   % For the 5 OVA classifiers
            
            % Corresponding cluster
            CL = cell2mat(Classif(1,C));
            REST = cell2mat(Classif(2,C));
            
            % Get feature values for Cluster and for rest of clusters
            FV_cluster = FeatVal(CL,:);
            FV_rest = FeatVal(REST,:);

            % Calculate sumrank for each feature
            for F = 1:numel(feat_id)                
                Pvalue(C,F) = ranksum(FV_cluster(:,F),FV_rest(:,F));           
            end
      
        end
        
        % Store in cell variable, for each iterations
        disp(sprintf('%s, %s iter',sub,string(Nf)))
        
        PVAL_CW(Nf,:) = Pvalue(1,:);
        PVAL_C1(Nf,:) = Pvalue(2,:);
        PVAL_C2(Nf,:) = Pvalue(3,:);
        PVAL_C3(Nf,:) = Pvalue(4,:);
        PVAL_CR(Nf,:) = Pvalue(5,:);

        
    end

   % Average across Nf iterations 
   PVAL_CW_iter(D,:) = mean(PVAL_CW);
   PVAL_C1_iter(D,:) = mean(PVAL_C1);
   PVAL_C2_iter(D,:) = mean(PVAL_C2);
   PVAL_C3_iter(D,:) = mean(PVAL_C3);
   PVAL_CR_iter(D,:) = mean(PVAL_CR);
   
   clear block
   
end

    % Average across D datasets
   PVAL_CW_iter_D = mean(PVAL_CW_iter);
   PVAL_C1_iter_D = mean(PVAL_C1_iter);
   PVAL_C2_iter_D = mean(PVAL_C2_iter);
   PVAL_C3_iter_D = mean(PVAL_C3_iter);
   PVAL_CR_iter_D = mean(PVAL_CR_iter);    
   

   % Re arrange into one array (Row 1 = C1, Row 5 = CR, Col = all F)
   PVAL_ALL_iter = [PVAL_CW_iter_D;PVAL_C1_iter_D;PVAL_C2_iter_D;PVAL_C3_iter_D;PVAL_CR_iter_D];
   
   
   % Sort features from best (lowest p) to worst (highest p) and get their
   % name / Operation type / pvalue
   
   for i = 1:5
       
       [~,I] = sort(PVAL_ALL_iter(i,:),'ascend');
       PVAL_ALL_iter(i,:) = PVAL_ALL_iter(i,I);
       
       % Get the 3 best features
       BestFeat = I(1:10);
       FeatName{i} = Operations.CodeString(BestFeat);
       FeatKey{i} = Operations.Keywords(BestFeat);
       
   end
       
       

 


%% No balanced dataset (Take all TS) (problem: using all epochs for training and testing: no cross validation 70/30, supervised

%%% --> OVA for each cluster, ranksum an ranking according to p value


%%% 1) Get the 5 clusters: (kmeans) and the feature values
%%% 2) Do sumrank for each feature, for each of the 5 OVA classifier
%%% 3) Rank fatures according to their log(p) value

Subs = {'439'}; % '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

%%%%% reduce to 5603 (for later)

load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Matrix_excl_all_feat_removed(5603)_all_datasets')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean_D_excl')
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

% Get 'Operations of WB features only, from HCTSA_N
load HCTSA.mat   
Operations_ID = [Operations.ID].';

equi_Top_Feat = setdiff(Operations_ID,spec_and_common_feat); % remove SV features
Operations = Operations(equi_Top_Feat,:);                      % 'Operations' with WB features only (from HCTSA_N)

CodeString = {Operations.CodeString}.';  
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

for D = 1:length(Subs)
    
    sub = Subs{D};

    % File selection
    WHICH_DATA = str2num(sprintf('%s',sub)); 

    % go to current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
   
    % Load TS_DataMat and Operations
    datam = load('HCTSA_N','TS_DataMat');
    load('HCTSA_N','Operations');
    all_op = Operations;

    %%%%% reduce to 5603 (get WB feature for that dataset)
    equi_Top_Feat = setdiff(Operations.ID,spec_and_common_feat); % remove SV features
    Idx_WB_Feat = find(ismember(Operations.ID, equi_Top_Feat));  % Index of WB features
    
    datamat = datam.TS_DataMat(:,Idx_WB_Feat);
    all_op = all_op(Idx_WB_Feat,:);

    % How many channels
    v = 1;   

    % Load data matrix for one feature
    feat_id = 1:5603;          % How many features (1)

    hctsa_ops = datamat(:,feat_id);       % hctsa_ops = datamat if # feature = 1

    %%%%%%%%%%%% Run k-means
    
    whichData = WHICH_DATA;
    annotation = load(sprintf('ccshs_1800%s_annot.mat',sub),'sleepstage');
    label = annotation.sleepstage;
    trainingProportion = 0.7;
    
    stgID = epochCounter(whichData,label);
    stgLab = {'W','N1','N2','N3','R'};

      % All TS from the dataset
     trainTS = 1:numel(label);
     testTS = 1:numel(label);

   
    if v==1
        trainMat = hctsa_ops(trainTS',:);
        testMat = hctsa_ops(testTS',:);
    else
        increment=round((size(hctsa_ops,1)/7))*(v-1);
        trainMat = [trainMat hctsa_ops((trainTS+increment)',:)];
        testMat = [testMat hctsa_ops((testTS+increment)',:)];
    end


     %%%%%%%%%%%%%%% k-means

     NUM_CLUSTERS = 5;

    [clustID,Kcentre,sse,Dist] = kmeans(trainMat,NUM_CLUSTERS,'Distance','sqeuclidean',...
            'Display','off','Replicates',50,'MaxIter',500);


     %%%%%%%%%%%%% Sequential maximum matching

     for n=1:length(testTS)

        % Calculate Euclidean distance from datapoint to each centre
        distance = sqrt(sum((bsxfun(@minus,Kcentre,testMat(n,:))).^2,2)); %%%% Calculate distance between each stage data point (Test mat) 
                                                                                    %%%% and center of cluster from kmeans. Data point with least distance from a center goes to cluster of that center. ClustTest = predictTest, = the cluster decisions
        % Find the cluster of minimum distance              
        [~,clustTest(n)] = min(distance);

     end


    %%% Find equivalent cluster-sleep stage pair

    stgLabel = [0,1,2,3,5];

    pro_no = zeros(length(stgID.nStg),length(unique(clustID)));

    for m = 1:length(unique(clustID))

        pro_id = find(clustID == m); 
        % Actual labelled stage
        actualStg = label(trainTS(pro_id));
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

    equi_stage = stgLabel(equi_class);   

    %%%%%%%%%%% Now, calculate ranksum

    %%%% 1) Get cluster names and the feature values

     CLUST = equi_stage(clustTest);   % Cluster mapped to stage

     % Name of cluster
     CLUST_NAME_idx = equi_class(clustTest);
     NAME = {'CW','C1','C2','C3','CR'};
     CLUST_NAME = NAME(CLUST_NAME_idx);

     % Feature values
     for i = 1:numel(testTS)
        FeatVal(i,:)= testMat(i,:);
     end


    %%%% 2) Calculate ranksum

    % Get index of epochs for each cluster
    CW = find(contains(CLUST_NAME,'CW'));
    C1 = find(contains(CLUST_NAME,'C1'));
    C2 = find(contains(CLUST_NAME,'C2'));
    C3 = find(contains(CLUST_NAME,'C3'));
    CR = find(contains(CLUST_NAME,'CR'));

    % Set the 5 OVA classifiers (IDs for feature values)
    CW_rest = [{CW};{[C1 C2 C3 CR]}];
    C1_rest = [{C1};{[CW C2 C3 CR]}];
    C2_rest = [{C2};{[C1 CW C3 CR]}];
    C3_rest = [{C3};{[C1 C2 CW CR]}];
    CR_rest = [{CR};{[C1 C2 C3 CW]}];

    Classif = [CW_rest C1_rest C2_rest C3_rest CR_rest];

    for C = 1:5   % For the 5 OVA classifiers

        % Corresponding cluster
        CL = cell2mat(Classif(1,C));
        REST = cell2mat(Classif(2,C));

        % Get feature values for Cluster and for rest of clusters
        FV_cluster = FeatVal(CL,:);
        FV_rest = FeatVal(REST,:);

        % Calculate sumrank for each feature
        for F = 1:numel(feat_id)                
            Pvalue(C,F) = ranksum(FV_cluster(:,F),FV_rest(:,F));           
        end

    end

    % Store in cell variable, for each iterations
    % fprintf('%s Iterations',string(Nf))
    PVAL_CW = Pvalue(1,:);
    PVAL_C1 = Pvalue(2,:);
    PVAL_C2 = Pvalue(3,:);
    PVAL_C3 = Pvalue(4,:);
    PVAL_CR = Pvalue(5,:);


   % Re arrange into one array (Row 1 = C1, Row 5 = CR, Col = all F)
   PVAL_ALL_iter = [PVAL_CW;PVAL_C1;PVAL_C2;PVAL_C3;PVAL_CR];

   % Sort features from best (lowest p) to worst (highest p) and get their
   % name / Operation type / pvalue

   for i = 1:5

       [~,I] = sort(PVAL_ALL_iter(i,:),'ascend');
       PVAL_ALL_iter(i,:) = PVAL_ALL_iter(i,I);

       % Get the 3 best features
       BestFeat = I(1:10);
       FeatName{i} = Operations.CodeString(BestFeat);
       FeatKey{i} = Operations.Keywords(BestFeat);

   end


end
 
      

      
