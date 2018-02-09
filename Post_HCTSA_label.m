%% Perform Post_HCTSA analysis with labeled Learn01
% Compare between supervised and unsupervised classification
homedir = pwd;
cd('.\081017_Learn02_Labeled')

% Load normalised HCTSA (Normalised seems to give better result)
load('HCTSA_N.mat') 
cd(homedir)

%% Re-order TS, such that each row represents featuresxchannels at each time segment
% delimiter =',';
Keywords = SUB_cell2cellcell({TimeSeries.Keywords}); % Split into sub-cells using comma delimiter

%% Obtain time label which is the second element in each cell
for keylab = 1:length(Keywords)
    % Find which component ID the name matched
    % Due to error when creating data matrix, after ICA the channel should
    % be called component instead, hence, channel F3 means component 1.
    timelabel(keylab) = Keywords{keylab,1}(2);
end

%% Group same time together
n=0;
timeuniq = unique(timelabel);
for time = 1:length(timeuniq) % 200:1250 % 1:length(timeuniq) %126:972 Learn02 %1:1090 Learn 01
    n=n+1;
    timename = sprintf('timeseg_%d',time);
    id = find(ismember(timelabel,timename));
    temp = TS_DataMat(id,:);
    % Reshape to make a row of the same time with all the features
    timedata(n,:) = reshape(temp,1,[]);
    % timetime{time} = timename;       
end
% Each column of timedata is the featurexchannel of time segments. Size of
% timedata is [number of time segment]x[number of features x channels]

%% Cluster time data using k-means clustering into 5 sleep stages + artifact stage
% Assume that each cluster correspond to each sleep stage
% Try different clustering algorithm BUT how to know which one is the best?
% 
nclust = 6; % Ideally, 5 sleep stages + 1 artifact

%% Additional step - detect wake vs sleep first and do clustering only for sleeping time
nWS = 2; 
clustWS = zeros(1,n);
[clustWS,centreWS,~] = kmeans(timedata,nWS,'Distance','sqeuclidean',...
                        'Display','final','Replicates',100,'MaxIter',1000);
%% Identify sleeping time
sleeptime = 99:969; % From previous clustering

%% K-means clustering
clustID = zeros(1,n);
[clustID(sleeptime),centrek,~] = kmeans(timedata(sleeptime,:),nclust,'Distance','sqeuclidean',...
                        'Display','final','Replicates',100,'MaxIter',500);

% FullClust(length(clustID))=struct; % Store time, clustering ID after processing and cause of NaN (in causeNaN struct)

%% Visualise clustering time course 
tt = 1:length(timeuniq); % Learn01 54:1074; Learn02 92:951; 1:length(timeuniq); % 200:1250;
figure;
stem(sleeptime,clustID(sleeptime),'Marker','.')
axis([0 n 0 nclust+1])
title('Clustering time course')
xlabel('Epoch (30-second segment)')
ylabel('Cluster #')

% %% tSNE - to visualise the clustering 
% addpath(genpath('F:/tSNE'));
% mapped_timedata = tsne(timedata,[],2,100,50);
% figure;
% gscatter(mapped_timedata(:,1),mapped_timedata(:,2),idk)
%% Transition probability matrix
% Obtain the transition probability matrix between clusters
% Hypothesis: Artifact clusters (Cluster 6) follows NREM > REM

% Call TPM function
% addpath(genpath('F:/AnalysisFunc/'))
TransMat = TPM(clustID(sleeptime),nclust);
TransGraph = digraph(TransMat);
figure;
plot(TransGraph,'EdgeLabel',TransGraph.Edges.Weight)
title('Transition Probability of clustering')
% To visualise the value of transition probability
figure;
colormap('gray')
imagesc(TransMat)
colorbar
%% Entropy calculation
% Probability of each stage P(X=x) where x is the cluster ID/stage
stageCount = zeros(1,6);
for t = sleeptime
    for i = 1:nclust
        if clustID(t)==i  % Count number of time spent on each stage
            stageCount(i)=stageCount(i)+1;
        end
    end
end
stageProb = stageCount./sum(stageCount);

% Calculate entropy
% H_dependent = entropy(TransMat,'conditional')
H_prior = -sum(stageProb.*log2(stageProb)); 

for i=1:nclust % For each current state, the probability of transitioning to next states is x(i,j)
    for j=1:nclust
        if TransMat(i,j)~=0    % If probability is non-zero, calculate p*log(p)
           H_transmat(i,j)= TransMat(i,j)*log2(TransMat(i,j));
           % H_ind(i,j)=x(i,j)*log2(x(i,i)/x(i,j));  % Value of each element
        else
           H_transmat(i,j)= 0;
        end
    end
    % Sum of probability in j (next stage | current stage)
   
end
% Entropy:
H_given = sum(H_transmat,2);
H_overall = -sum(stageProb.*H_given');
TPMinfo = H_prior - H_overall;
%% Save cluster ID to compare
% cd('./020917_Learn03')
% save('270917_Learn03_Result','tt','clustID','TransMat',...
%      'stageProb','stageCount','H_prior','H_overall','TPMinfo')
