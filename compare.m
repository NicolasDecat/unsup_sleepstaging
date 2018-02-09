% Compare scored data and analysis result data
%% Load data
homedir = pwd;
cd('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\310817\learn\polysomnography\annotations-events-profusion')
load('Learn03Annot.mat')
cd('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\020917_Learn03')
load('270917_Learn03_Result.mat')
cd(homedir)

%% Transpose sleepstage data
% Same size
sleepstage = sleepstage';

%% Re-assign clustering
% To compare the stages, align sleep part (omit wake stage before and
% after)
sleeptime = 112:912; % Learn01 54:1074; Learn02 92:951;% 112:912;
scored = sleepstage(sleeptime);
result = clustID(sleeptime);
nclust = 6;
%% Using probability to compare
% Probability of each stage P(X=x) where x is the cluster ID/stage
stageCount = zeros(2,6);
for t = sleeptime
    for i = 1:nclust
        if clustID(t)==i  % Count number of time spent on each stage
            stageCount(1,i)=stageCount(1,i)+1;
        end
        if sleepstage(t)==(i-1)
            stageCount(2,i)=stageCount(2,i)+1;
        end
    end
end
stageProb = stageCount./sum(stageCount,2);

%% which one is similar? 
[sProb,stageID]=sort(stageProb,2);

%%  Visualise
n=length(sleepstage);
figure;
subplot(2,1,1)
plot(1:n,clustID,'r')
axis([0 1200 0 6])
legend('Result')
subplot(2,1,2)
plot(1:n,sleepstage)
axis([0 1200 0 6])
legend('Scored')
%% Change stage in result data to match scored data
% For comparison
for t=1:length(result)
    matchID = find(stageID(1,:)==result(t));
    resultN(t)=stageID(2,matchID)-1; % Back to 0-5 stages
end

figure;
subplot(3,1,1)
plot(1:n,sleepstage)
axis([0 1200 0 6])
legend('scored')
subplot(3,1,2)
plot(sleeptime,resultN)
axis([0 1200 0 6])
legend('result')

%% Calculate percentage accuracy
correct = 0; wrong = 0;
for i=1:length(resultN)
    if resultN(i)==scored(i)
        correct = correct+1;
        resultID(i) = 1;
    else
        wrong = wrong+1;
        resultID(i) = NaN;
    end
end
correctPer = correct/(correct+wrong)
subplot(3,1,3)
plot(sleeptime,resultN)
legend('result with correct points')
hold on
correctpoints = resultID.*resultN;
plot(sleeptime,correctpoints,'rx')
axis([0 1200 0 6])
hold off

%% Correlation
clear r
r = corrcoef(resultN,scored);

figure;
imagesc(r)
colormap gray
colorbar