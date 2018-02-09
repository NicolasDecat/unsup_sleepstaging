% Compare scored data and analysis result data
%% Load data
homedir = pwd;
cd('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\310817\learn\polysomnography\annotations-events-profusion')
load('Learn02Annot.mat')
% cd('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\020917_Learn03')
% load('learn01_clustering_WS.mat')
cd(homedir)

%% Transpose sleepstage data
% Same size
sleepstage = sleepstage';

%% Re-assign clustering
% To compare the stages, align sleep part (omit wake stage before and
% after)
% sleeptime = 1:length(sleepstage); % Learn01 54:1074; Learn02 92:951;% Learn03: 112:912;
scored = sleepstage(sleeptime);
result = clustID(sleeptime);
nclust = 6;
%% Using probability to compare
% Probability of each stage P(X=x) where x is the cluster ID/stage
stageCount = zeros(2,6);
for t = 1:length(sleeptime)
    for i = 1:nclust
        if result(t)==i  % Count number of time spent on each stage
            stageCount(1,i)=stageCount(1,i)+1;
        end
        if scored(t)==(i-1)
            stageCount(2,i)=stageCount(2,i)+1;
        end
    end
end
stageProb = stageCount./sum(stageCount,2);

%% which one is similar? 
[sProb,stageID]=sort(stageProb,2);

%%  Visualise
% n=length(sleepstage);
figure;
subplot(2,1,1)
plot(sleeptime,result,'r')
axis([0 n 0 6])
legend('Result')
subplot(2,1,2)
plot(1:n,sleepstage)
axis([0 n 0 6])
legend('Scored')
%% Change stage in result data to match scored data
% For comparison
for t=1:length(result)
    matchID = find(stageID(1,:)==result(t));
    resultN(t)=stageID(2,matchID)-1; % Back to 0-5 stages
end
%%
camelot = [0.5117,0.2109,0.2656];
adjCamelot = [160/255,73/255,77/255];
genoa = [0.2109,0.5117,0.4570];
adjGenoa = [53/255,155/255,126/255];
eucalyptus = [0.2031,0.6016,0.4023];

% Only for learn02
for s=1:length(resultN)
    if resultN(s)==4 % If stage is N4, change it to N3
        resultN5(s)=3;
    else
        resultN5(s)=resultN(s);
    end
end
% ============
resultNN = zeros(1,length(sleepstage));
resultNN(sleeptime)=resultN5;
figure;
subplot(3,1,1)
plot(1:n,sleepstage,'Color',adjCamelot)
axis([0 length(sleepstage) 0 6])
legend('scored')
subplot(3,1,2)
plot(1:n,resultNN,'Color',adjGenoa)
axis([0 n 0 6])
legend('result')

%% Calculate percentage accuracy
correct = 0; wrong = 0;

for i=1:length(resultNN)
    if resultNN(i)==sleepstage(i)
        correct = correct+1;
        resultID(i) = 1;
    else
        wrong = wrong+1;
        resultID(i) = NaN;
    end
end
correctPer = correct/(correct+wrong)
subplot(3,1,3)
plot(1:n,resultNN,'Color',adjGenoa)
legend('result with correct points')
hold on
correctpoints = resultID.*resultNN;
plot(1:n,correctpoints,'Color',adjCamelot,'Marker','d','MarkerSize',4)
axis([0 n 0 6])
hold off

%% Correlation
r_aligned = corrcoef(resultN,scored);
corr_aligned = r_aligned(2,1);
r_raw = corrcoef(result,scored);
corr_raw = r_raw(2,1);

%% Classification confusion matrix
addpath(genpath('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging\HCTSA'));

% Modified from TS_classify (Ben's)
scoredLabel = BF_ToBinaryClass(sleepstage,n_clust);
predictLabel = BF_ToBinaryClass(resultNN,n_clust);
figure
plotconfusion(scoredLabel,predictLabel);

% Fix axis labels:
n_clust = 5;
ax = gca;
ax.XTickLabel(1:n_clust) = groupNames;
ax.YTickLabel(1:n_clust) = groupNames;
ax.TickLabelInterpreter = 'none';

% Make a nice white figure background:
f = gcf; f.Color = 'w';



%% Confusion matrix (modified)
% [C,order] = confusionmat(scored,resultN)
% %%
% n = 4;
% m = 5;
% Amax = 7;
% cmap = autumn(Amax);
% A = randi(Amax,m,n);
% figure
% image(A)
% colormap(cmap)
% set(gca,'XTick',[],'YTick',[],'YDir','normal')
% [x,y] = meshgrid(1:n,1:m);
% text(x(:),y(:),num2str(A(:)),'HorizontalAlignment','center')
% 
% %%
% M=[1,2,3;4,3,2;1,4,2]
% figure
% imagesc(M)
% for a=1:3
%     for b=1:3
%         text(b,a,num2str(M(a,b)))
%     end
% end