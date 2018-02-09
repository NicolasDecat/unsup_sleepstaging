%% Cross-validation code
% Example: learn01 data
% After reading the annotation file from xml using read_annot.m
%addpath(genpath('/Users/sleeping/Documents/MATLAB/ccshs_data'))
annotation = load('ccshs_1800001_annot.mat');
label = annotation.sleepstage;

%% Proportion of each sleep stage (0 - wake, 1-4 NREM, 5 - REM)
stgNum = size(unique(label));
stgLab = {'W','N1','N2','N3','N4','R'}; % {'W','N1','N2','N3','R'};

%======== Just to check but needed? =========
for i = 1:stgNum
    stgID.stgPro(i) = sum(label==i-1);
end
% ===============

% Sleep stage IDs
w = 0; n1 = 0; n2 = 0; n3 = 0; n4 = 0; r = 0;
for n = 1:length(label)
    switch label(n)
        case 0 % Wake
            w=w+1;
            stgID.allID.W(w) = n;
        case 1 % N1
            n1=n1+1;
            stgID.allID.N1(n1) = n;
        case 2 % N2
            n2=n2+1;
            stgID.allID.N2(n2) = n;
        case 3 % N3
            n3=n3+1;
            stgID.allID.N3(n3) = n;
        case 4 % N4
            n4=n4+1;
            stgID.allID.N4(n4) = n;
        case 5 % REM
            r=r+1;
            stgID.allID.R(r) = n;
    end
end

stgID.stgPro = [w, n1, n2, n3, n4, r];
% ===============

clear w n1 n2 n3 n4 r
%% Remove class with less than cut-off
cutoff = 0.02*length(label);
n=0;
for i = 1:length(stgID.stgPro)
    if stgID.stgPro(i)>=cutoff
        n=n+1;
        stgID.useStg(n)=i;
        stgID.usePro(n)=stgID.stgPro(i);
        stgID.useStgName(n) = stgLab(i);
    end
end
%% Minimum samples
stgID.Nmin = min(stgID.usePro);

%% Multiple iteration of randomisation and cross-validation
for Nf = 1:100
%% Random sampling from each class
train70 = round(0.7*stgID.Nmin);

for m=1:length(stgID.useStg)
    randID = randperm(stgID.usePro(m),stgID.Nmin);
    allID = stgID.allID.(stgLab{stgID.useStg(m)});
    useID = allID(randID);
    stgID.useID.(stgLab{stgID.useStg(m)}) = useID;
    % 70% training
    randtrain = randperm(stgID.Nmin);
    trainID = useID(randtrain(1:train70));
    testID = useID(randtrain((train70+1):end));
    stgID.trainID.(stgLab{stgID.useStg(m)})= trainID;
    stgID.testID.(stgLab{stgID.useStg(m)})= testID;
    % Combine training ID to perform classification
    if m==1 % First iteration
        trainTS = trainID;
        testTS = testID;
    else
        trainTS = [trainTS,trainID]; 
        testTS = [testTS,testID];
    end
    clear randID allID useID randtrain trainID testID
end
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

% %% tSNE
% % Set parameters
% no_dims = 2; 
% initial_dim = 50;
% perplexity = 50;
% 
% % Run t-SNE
% mapped_data = tsne(trainMat,[],no_dims,initial_dim,perplexity); % Unlabelled data
% 
% % Plot result
% figure;
% subplot(1,2,1)
% gscatter(mapped_data(:,1),mapped_data(:,2),label(trainTS));
% title('Labelled')
% subplot(1,2,2)
% gscatter(mapped_data(:,1),mapped_data(:,2),clustID)
% title('Clustered')
% %% Visualise labelled and clustered of training data
% [sortID, I] = sort(testTS);
% figure;
% plot(sortID,clustTest(I),sortID,label(testTS(I)))
% %%
% figure;
% plot(trainTS,clustID,'o',testTS,label(testTS),'*')
% hold on
% plot(1:1374,label)
% hold off

%% Count number of sample from each class that are assigned to THE prototype (cluster)
% Highest count = equivalent cluster
for m = 1:length(unique(clustID))
    % ID that belongs to this cluster
    pro_id = find(clustID==m);
    N = histcounts(label(trainTS(pro_id)))
    
    for n = 1:stgNum
        pro_no(n) = sum(label(trainTS(pro_id))==stgID.useStg(n)-1);
    end
    disp(pro_no)
    [~,equi_class(m)]=max(pro_no);
    % Change to equivalent stage ID (0,1,2,3,5)
    
end
equi_stage = stgID.useStg(equi_class)-1;

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

end % End Nf-th randomisation

%% Average output
Output(k).testCorrect = mean(Pcorrect);
Output(k).testError = mean(Perror);
Output(k).trainCorrect = mean(trainPcorrect);
Output(k).trainError = mean(trainPerror);


%% Confusion matrix - for both train and test
% For one run(k)... can change which data to use later
addpath(genpath('/Users/sleeping/Documents/MATLAB/unsup_sleep_staging/HCTSA'))


%% Confusion matrix of train data

% Labelled - make non-zero stage
g_labelTrain = label(trainTS)+1;

% Clustered - Use final clustering output
g_clustTrain = equi_train+1;

% Cluster 6 becomes 5
g_labelTrain(g_labelTrain==6) = 5;
g_clustTrain(g_clustTrain==6) = 5;

% BINARY TO CLASS FUNCTION FROM BEN'S HCTSA
labelTrainBF= BF_ToBinaryClass(g_labelTrain,nclust);
clustTrainBF = BF_ToBinaryClass(g_clustTrain,nclust);


% Visualise confusion matrix
figure;
plotconfusion(labelTrainBF,clustTrainBF)

% Plot setting
ax = gca;
ax.XTickLabel(1:nclust)=stgID.useStgName;
ax.YTickLabel(1:nclust)=stgID.useStgName;


%% Confusion matrix of test data
% Labelled - make non-zero stage
g_labelTest = label(testTS)+1;

% Clustered - Use final clustering output
g_clustTest = equi_test+1;

% Cluster 6 becomes 5
g_labelTest(g_labelTest==6) = 5;
g_clustTest(g_clustTest==6) = 5;

% BINARY TO CLASS FUNCTION FROM BEN'S HCTSA
labelTestBF= BF_ToBinaryClass(g_labelTest,nclust);
clustTestBF = BF_ToBinaryClass(g_clustTest,nclust);


% Visualise confusion matrix
figure;
plotconfusion(labelTestBF,clustTestBF)

% Plot setting
ax = gca;
ax.XTickLabel(1:nclust)=stgID.useStgName;
ax.YTickLabel(1:nclust)=stgID.useStgName;


%% Clear variables for the next run
clearvars -except Output datamat feat_id features k complexity