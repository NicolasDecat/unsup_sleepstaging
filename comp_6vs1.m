%% Compare clustering of 6 channels and 1 channel 
% Find the features which describe the sleep stage
homedir = pwd;
% Load data from 1 channel
cd('./200717')
data_one = load('postHCTSA1chan_topops.mat');
cd(homedir)
clustID_one = data_one.clustID;
% Load data from 6 channels
cd('./270817')
data_six = load('postHCTSA6chan_topops.mat');
cd(homedir)
clustID_six = data_six.clustID;
%% Plot clustering to compare
figure;
ax1 = subplot(2,1,1);
plot(clustID_one)
axis([0 length(clustID_one) 0 7])
legend('1 channel (Fz)')
ax2 = subplot(2,1,2);
plot(clustID_six)
axis([0 length(clustID_six) 0 7])
legend('6 channels')
linkaxes([ax1,ax2],'x')
%% Compare two datasets
sleeptime = 1:length(clustID_one);
nclust = data_one.nclust;
% Probability of each stage (P(X=x) where x is the cluster ID/stage
stageCount = zeros(2,6);
for t = sleeptime
    for i = 1:nclust
        if clustID_one(t)==i  % Count number of time spent on each stage
            stageCount(1,i)=stageCount(1,i)+1;
        end
        if clustID_six(t)==i    % 2nd row = 6 channel data
            stageCount(2,i)=stageCount(2,i)+1;
        end
    end
end
stageProb = stageCount./sum(stageCount,2);
%% Align stage according to probability
[sProb,stageID]=sort(stageProb,2);
% For comparison
for t=1:length(clustID_one)
    matchID = find(stageID(1,:)==clustID_one(t));
    clustID_oneN(t)=stageID(2,matchID); % Back to 0-5 stages
end
%% Plot clustering to compare
figure;
ax1 = subplot(3,1,1);
plot(clustID_oneN)
axis([0 length(clustID_one) 0 7])
legend('1 channel (Fz) - aligned')
ax2 = subplot(3,1,2);
plot(clustID_six)
axis([0 length(clustID_six) 0 7])
legend('6 channels')
linkaxes([ax1,ax2],'x')
%% Calculate percentage accuracy
correct = 0; wrong = 0;
for i=1:length(clustID_oneN)
    if clustID_oneN(i)==clustID_six(i)
        correct = correct+1;
        resultID(i) = 1;
    else
        wrong = wrong+1;
        resultID(i) = NaN;
    end
end
correctPer = correct/(correct+wrong)
ax3=subplot(3,1,3);
plot(sleeptime,clustID_oneN)
legend('result with correct points')
hold on
correctpoints = resultID.*clustID_oneN;
plot(sleeptime,correctpoints,'rx')
axis([0 length(clustID_oneN) 0 7])
hold off
%% Correlation
r = corrcoef(clustID_oneN,clustID_six);
correlation = r(2,1)
%% Features
% Compare accuracy of each features in 1channe dataset with 6 channel
TS_one = data_one.TS_DataMat;
center_one = data_one.centrek;
TS_six = data_six.TS_DataMat;
TS_sixR = data_six.timedata;
%% Clustering using one feature - Defining feature
oneChan = struct;

for f=1:195 % Go through each feature
    % Compare value of time feature at each time segment with the centre of
    % cluster
    oneChan(f).timexfeat = TS_one(:,f);
    oneChan(f).centrefeat = center_one(:,f);
    for t=1:length(TS_one(:,f))
        % Distance from center of each cluster
        oneChan(f).centre_dist(:,t) = abs(oneChan(f).centrefeat - ((oneChan(f).timexfeat(t))*ones(6,1)));
        % Closest cluster
        [~,clustxfeat(t,f)]=min(oneChan(f).centre_dist);
        % Closest cluster -> align with 6 chan data
        matchID = find(stageID(1,:)==clustxfeat(t,f));
        clustxfeatAlign(t,f)=stageID(2,matchID); % Back to 0-5 stages
        
        % Count number of stages -> Probability of sleep stage
        stageCountxfeat=zeros(1,6);
        for i=1:nclust
            if clustxfeatAlign(t,f)==i
                stageCountxfeat(i)=stageCountxfeat(i)+1;
            end
        end
    end
     oneChan(f).clustxfeat=clustxfeatAlign(:,f );
    % Compare with 6 channel data
    % Entropy calculation
    oneChan(f).stageProbxfeat = stageCountxfeat./sum(stageCountxfeat);
    oneChan(f).TransMat = TPM(clustxfeatAlign(:,f),nclust);
    oneChan(f).H_prior = -sum(oneChan(f).stageProbxfeat.*log2(oneChan(f).stageProbxfeat)); 
    for i=1:size(oneChan(f).TransMat,1) % For each current state, the probability of transitioning to next states is x(i,j)
        for j=1:size(oneChan(f).TransMat,2)
            if oneChan(f).TransMat(i,j)~=0    % If probability is non-zero, calculate p*log(p)
                H_transmat(i,j)= oneChan(f).TransMat(i,j)*log2(oneChan(f).TransMat(i,j));
           % H_ind(i,j)=x(i,j)*log2(x(i,i)/x(i,j));  % Value of each element
            else
                H_transmat(i,j)= 0;
            end
        end
    end
   % Entropy:
   H_given = sum(H_transmat,2);
   oneChan(f).H_overall = -sum(oneChan(f).stageProbxfeat.*H_given');
   oneChan(f).TPMinfo = oneChan(f).H_prior - oneChan(f).H_overall;
    % Calculate percentage accuracy - assuming 6channel data is the most
    % accurate
    correct = 0; wrong = 0;
    for t=1:length(TS_one(:,f))
        if clustID_oneN(i)==clustID_six(i)
            correct = correct+1;
            oneChan(f).resultID(i) = 1;
        else
            wrong = wrong+1;
            oneChan(f).resultID(i) = NaN;
        end
    end
    oneChan(f).correctTimeSeg = correct;
    oneChan(f).wrongTimeSef = wrong;
    oneChan(f).correctPer = correct/(correct+wrong);
    % Correlation
    rho = corrcoef(clustxfeatAlign(:,f),clustID_six);
    oneChan(f).correlation = rho(2,1);
    fprintf('Finish feature #%d\n',f)
end
% Combinations of feature (later)