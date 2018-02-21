% %% Cross-validation code
% configuration_settings; 
% 
% %% K-means clustering + Nearest centroid classifier
% % Example: learn01 data
% % After reading the annotation file from xml using read_annot.m
% %addpath(genpath('/Users/sleeping/Documents/MATLAB/ccshs_data'))
% annotation = load(ANSWER_FILE);
% label = annotation.sleepstage;
% 
% %% Proportion of each sleep stage (0 - wake, 1-4 NREM, 5 - REM)
% stgNum = size(unique(label));
% stgLab = {'W','N1','N2','N3','N4','R'}; % {'W','N1','N2','N3','R'};
% 
% %% Remove initial W stage from randomisation and sampling
% % Marking the end of W stage
% endW = [334,380,391,375,174];
% endS = [1374,1442,1442,1531,1492];
% 
% % endID = find(whichData==validData);
% endID = 1;
% whichData = 1;
% if whichData == 14 % Remove the ending epochs for the ccshs_1800014
%     selectID = [(endW(endID)+1):end14(1)]';
% else
%     selectID = [endW(endID)+1:length(label)]';
% end
% selectLabel = label(selectID);
% 
% %% Counting number of each stage
% 
% %======== Just to check but needed? =========
% for i = 1:stgNum
%     stgID.stgPro(i) = sum(label==i-1);
% end
% % ===============
% 
% % Sleep stage IDs
% w = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0; r = 0;
% for n = 1:length(selectLabel)
%     switch selectLabel(n)
%         case 0 % Wake
%             w=w+1;
%             stgID.allID.W(w) = n;
%         case 1 % N1
%             n1=n1+1;
%             stgID.allID.N1(n1) = n;
%         case 2 % N2
%             n2=n2+1;
%             stgID.allID.N2(n2) = n;
%         case 3 % N3
%             n3=n3+1;
%             stgID.allID.N3(n3) = n;
%         case 4 % N4
%             n4=n4+1;
%             stgID.allID.N4(n4) = n;
%         case 5 % REM
%             r=r+1;
%             stgID.allID.R(r) = n;
%     end
% end
% 
% stgID.stgPro = [w, n1, n2, n3, n4, r];
% % ===============
% 
% clear w n1 n2 n3 n4 r
% %% Remove class with less than cut-off
% cutoff = 0.02*length(label);
% n=0;
% for i = 1:length(stgID.stgPro)
%     if stgID.stgPro(i)>=cutoff
%         n=n+1;
%         stgID.useStg(n)=i;
%         stgID.usePro(n)=stgID.stgPro(i);
%         stgID.useStgName(n) = stgLab(i);
%     end
% end
% %% Minimum samples
% stgID.Nmin = min(stgID.usePro);

%% Multiple iteration of randomisation and cross-validation
for Nf = 1:20
%% Random sampling from each class
% train70 = round(0.7*stgID.Nmin);
% 
% for m=1:length(stgID.useStg)
%     randID = randperm(stgID.usePro(m),stgID.Nmin);
%     allID = stgID.allID.(stgLab{stgID.useStg(m)});
%     useID = allID(randID);
%     stgID.useID.(stgLab{stgID.useStg(m)}) = useID;
%     % 70% training
%     randtrain = randperm(stgID.Nmin);
%     trainID = useID(randtrain(1:train70));
%     testID = useID(randtrain((train70+1):end));
%     stgID.trainID.(stgLab{stgID.useStg(m)})= trainID;
%     stgID.testID.(stgLab{stgID.useStg(m)})= testID;
%     % Combine training ID to perform classification
%     if m==1 % First iteration
%         trainTS = selectID(trainID);
%         testTS = selectID(testID);
%     else
%         trainTS = [trainTS;selectID(trainID)]; 
%         testTS = [testTS;selectID(testID)];
%     end
%     clear randID allID useID randtrain trainID testID
% end
% trainTS and testTS are the time segment ID for training set and test set
% respectively.

%% Select data of wanted time ID
% **** Use timedata for >1 channel data. For HCTSA_N.mat of 1 channel, the
% output is the same.

%% Clustering with training dataset
% featrange = : ; % For now
trainMat = hctsa_ops(trainTS,:); % Change the number of features used*****
testMat = hctsa_ops(testTS,:);
nclust =5;
[clustID,centre,~] = kmeans(trainMat,nclust,'Distance','sqeuclidean',...
                        'Display','final','Replicates',50,'MaxIter',500);
                    
%% Centre of sleep stages
% Mean feature value of the labelled time segment
for m = 1:length(stgID.useStg)
    stgTrainID = stgID.trainID.(stgLab{stgID.useStg(m)});
    stgTrain = hctsa_ops(stgTrainID,:);
    stgCntr(m,:) = mean(stgTrain);
end
% Centre of clusters = centre

%% Test set classification (Nearest centroid classifier)
% using centre/mean features - minimum distance
for n=1:length(testTS)
    distance = sqrt(sum((centre-testMat(n,:)).^2,2));
    distance_stg = sqrt(sum((stgCntr-testMat(n,:)).^2,2));
    
%     distance = sum(sqrt((centre-repmat(testMat(n,:),4,1)).^2),2);
    [~,clustTest(n)] = min(distance);
    [~,clustTest_stg(n)] = min(distance_stg);
end
clustTest = clustTest';

%% tSNE
% 
% addpath(genpath('/Volumes/Seagate Expansion Drive/tSNE/'));
% % Set parameters
% no_dims = 2; 
% initial_dim = 50;
% perplexity = 50;
% 
% % Perform t-SNE for training and test together
% 
% % Concatenate training and test data matrix
% concatMat = [trainMat;testMat];
% 
% 
% % Run t-SNE
% mapped_data = tsne(concatMat,[],no_dims,initial_dim,perplexity); % Unlabelled data
% 
% % Plot result
% figure;
% subplot(1,2,1)
% gscatter(mapped_data(1:130,1),mapped_data(1:130,2),label(trainTS),[],'.....',16)
% hold on
% gscatter(mapped_data(131:185,1),mapped_data(131:185,2),label(testTS),[],'xxxxx',8)
% legend('W','N1','N2','N3','R')
% %h=plot(NaN,NaN,'k.',NaN,NaN,'kx')
% %legend(h,'Training','Test')
% hold off
% title('Labelled')
% axis([-10 6 -8 8])
% 
% subplot(1,2,2)
% gscatter(mapped_data(131:185,1),mapped_data(131:185,2),clustTest,[],'xxxxx',8)
% hold on
% gscatter(mapped_data(1:130,1),mapped_data(1:130,2),clustID,[],'.....',16)
% hold off
% title('Clustered')
% axis([-10 6 -8 8])
%% Visualise labelled and clustered of training data
[sortID, I] = sort(testTS);
figure;
plot(sortID,clustTest(I),sortID,label(testTS(I)))
%%
figure;
plot(trainTS,clustID,'o',testTS,label(testTS),'*')
hold on
plot(1:1374,label)
hold off

%% Count number of sample from each class that are assigned to THE prototype (cluster)
% Highest count = equivalent cluster
for m = 1:length(unique(clustID))
    % ID that belongs to this cluster
    pro_id = find(clustID==m);
    N = histcounts(label(trainTS(pro_id)));
    
    for n = 1:stgNum
        pro_no(m,n) = sum(label(trainTS(pro_id))==stgID.useStg(n)-1);
    end
    disp(pro_no);
    [~,equi_class(m)]=max(pro_no(m,:));
    % Change to equivalent stage ID (0,1,2,3,5)
    
end
equi_stage = stgID.useStg(equi_class)-1;

%% Equivalent stage - must be unique *** FIX THIS!!!

% if length(unique(equi_stage))~=length(unique(clustID))
%     % If some stages repeat, the one with more number of members would be
%     % assigned as that stage. The other stage would be the remaining
%     % stage.
%     N = histcounts(equi_stage); % Count the number of each stage, if all are uniqe, count = 1.
%     repeatIndex = find(N~=1);
%     repeatStg = N(repeatIndex(1));
%     % Looking at pro_no to determine which stage corresponds to the cluster
%     pro_repeat = pro_no(repeatstage,:);
%     % The one with most number of ... gets the stage
%     
% end
    
%% Equivalent stage - k-NN means
for a = 1:length(clustTest_stg)
    if clustTest_stg(a)==5
        clustTest_stg(a)=5;
    else
        clustTest_stg(a)=clustTest_stg(a)-1;
    end
end

correct_stg = sum(clustTest_stg==label(testTS)')

%% Training error
% Convert clustID into equivalent sleep stage label
for i=1:length(clustID)
    equi_train(i) = equi_stage(clustID(i));
end

%% Test error
% Convert clustTest into equivalent stage label
for i=1:length(clustTest)
    equi_test(i) = equi_stage(clustTest(i));
end

%% Percentage correct
trainCorrect = sum(equi_train==label(trainTS)');
trainIncorrect = sum(equi_train~=label(trainTS)');
testCorrect = sum(equi_test==label(testTS)');
testIncorrect = sum(equi_test~=label(testTS)');

% correct = sum(equi_test==label(testTS)');
% incorrect = sum(equi_test~=label(testTS)');
Pcorrect(Nf) = testCorrect/(testCorrect+testIncorrect);
Perror(Nf) = testIncorrect/(testCorrect+testIncorrect);
trainPcorrect(Nf) = trainCorrect/(trainCorrect+trainIncorrect);
trainPerror(Nf) = trainIncorrect/(trainCorrect+trainIncorrect);

%% Confusion matrix input
scoredTrain(Nf,:) = label(trainTS);
scoredTest(Nf,:) = label(testTS);

predictTrain(Nf,:)= equi_train;
predictTest(Nf,:)=equi_test;


end % End Nf-th randomisation

%% Average output
Output(k).testCorrect = mean(Pcorrect);
Output(k).testError = mean(Perror);
Output(k).trainCorrect = mean(trainPcorrect);
Output(k).trainError = mean(trainPerror);
%% Hypnogram
figure;
plot(trainTS,equi_train,'o',testTS,equi_test,'*')
hold on
plot(1:1374,label)
hold off
legend('Training','Test')

%% Confusion matrix - for both train and test
% For one run(k)... can change which data to use later
addpath('/Users/sleeping/Documents/MATLAB/unsup_sleep_staging/HCTSA/PeripheryFunctions/BF_ToBinaryClass.m')


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

% BINARY TO CLASS FUNCTION FROM BEN'S HCTSA
labelTrainBF= BF_ToBinaryClass(g_labelTrain,nclust);
clustTrainBF = BF_ToBinaryClass(g_clustTrain,nclust);



% Visualise confusion matrix
figure;
plotconfusion_custom(g_labelTrain, g_clustTrain, 'Confusion Matrix - Training');
saveas(gcf, strcat(CM_SAVE_DIR, filesep, 'CM_TRN_', int2str(k), '.png'));
%plotconfusion(labelTrainBF,clustTrainBF)

% Plot setting
% ax = gca;
% ax.XTickLabel(1:nclust)=stgID.useStgName;
% ax.YTickLabel(1:nclust)=stgID.useStgName;


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

% BINARY TO CLASS FUNCTION FROM BEN'S HCTSA
labelTestBF= BF_ToBinaryClass(g_labelTest,nclust);
clustTestBF = BF_ToBinaryClass(g_clustTest,nclust);


% Visualise confusion matrix
plotconfusion_custom(g_labelTest, g_clustTest, 'Confusion Matrix - Testing');
saveas(gcf, strcat(CM_SAVE_DIR, filesep, 'CM_TST_', int2str(k), '.png'));

% Plot setting
% ax = gca;
% ax.XTickLabel(1:nclust)=stgID.useStgName;
% ax.YTickLabel(1:nclust)=stgID.useStgName;


%% Clear variables for the next run
clearvars -except Output datamat feat_id features k complexity