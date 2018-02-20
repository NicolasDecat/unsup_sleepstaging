% Post HCTSA
% Following HCTSA analysis and combineBatchFiles.m 
configuration_settings

%% Normalise HCTSA.m
homedir = pwd;
% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)
% % Load data and normalise
cd(HCTSA_DATA_DIR)
TS_normalize('scaledRobustSigmoid',[0.8,1.0]);
% % TS_LabelGroups({'T1','T2','T3','T4'},'raw');
TS_plot_DataMatrix('norm');

%% Load normalised data matrix
load('HCTSA_N.mat')
cd(homedir)
% Re-order TS, such that each row represents featuresxchannels at each time
% segment
delimiter =',';
Keywords = SUB_cell2cellcell({TimeSeries.Keywords},delimiter); % Split into sub-cells using comma delimiter

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
%% K-means clustering
clustID = zeros(1,n);
[clustID,centrek,~] = kmeans(timedata,nclust,'Distance','sqeuclidean',...
                        'Display','final','Replicates',50,'MaxIter',500);

% FullClust(length(clustID))=struct; % Store time, clustering ID after processing and cause of NaN (in causeNaN struct)

%% Visualise clustering time course 
tt = 1:length(timeuniq); % Learn01 54:1074; Learn02 92:951; 1:length(timeuniq); % 200:1250;
figure;
stem(tt,clustID,'Marker','.')
axis([0 length(timeuniq) 0 nclust+1])
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
TransMat = TPM(clustID,nclust);
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
for t = tt
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
%     'stageProb','stageCount','H_prior','H_overall','TPMinfo')
 
 % %% Remove time segment with artifacts (cluster #6)**** (Arbitrary - 
% % need to find some better way of finding artifact cluster)
% % Cluster 6 seems to be cluster of time segment with large, visible
% % artifacts (Large motion, abnormally high amplitude)
% % By removing these cluster from the whole timecourse clustering plot, it
% % might reveal better clustering dynamic
% 
% % Create a struct that will keep tract of cause of NaN(why NaN is assigned
% % to this particular time segment)
% 
% for time= 1:length(clustID)
%     % Unprocessed fields are time and Orgcluster
%     FullClust(time).Timeseg = time;
%     FullClust(time).OrgCluster = clustID(time);
%     
%     if (FullClust(time).OrgCluster == 6) % If the it is artifact cluster (use struct for consistency)
%         FullClust(time).FltCluster = NaN;  % filtered Cluster ID
%         FullClust(time).CauseNaN.Artifact = 1;
%         % Replaced with NaN or adjacent cluster? 
%     else
%         FullClust(time).FltCluster = FullClust(time).OrgCluster ;
%         FullClust(time).CauseNaN.Artifact = 0;
%     end
% end
% 
% figure
% stem([FullClust.Timeseg],[FullClust.FltCluster],'Marker','.')
% axis([0 790 0 nclust+1])
% title('Clustering time course (with artifact filtering)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')

%%  3. Predict the clustering of candidate time segment using adjacent clustering information
% By looking at FltCluster to predict the cluster of candidate time
% segments
% 
% Fullclust.Predicted
% for time=1:length(clustID)
%     Extra condition for first and last time segments - indecisive time
%     segment
%     if (time==1)||(time==length(clustID))
%         FullClust(time).Predicted = NaN;
%         FullClust(time).CauseNaN.Edge= 1;
%         continue;
%     end
%     Conditions/ Rules for assigning cluster ID
%     1.) If the two adjacent time segment have the same cluster ID, then
%     the middle one also have the same cluster
%     if FullClust(time-1).FltCluster == FullClust(time+1).FltCluster % If time is between time of stable stages (not during transition) 
%        FullClust(time).Predicted = FullClust(time-1).FltCluster;
%        FullClust(time).CauseNaN.Edge = 0;      % Transition between odd time
%     2.) If the cluster ID of previous time NaN but the cluster of
%     next time exists,...  
%     elseif (isnan(FullClust(time-1).FltCluster))&&(~isnan(FullClust(time+1).FltCluster)) %If the previous time segment is NaN, then look at the one before it
%         Look at the time segment before that to see if they are the
%         same.
%         if FullClust(time+1).FltCluster == FullClust(time-2).FltCluster
%            FullClust(time).Predicted = FullClust(time+1).FltCluster; % time+1 cluster is assigned in previous iteration
%            FullClust(time).CauseNaN.Edge = 0;
%         else
%             If the time segment is at the transition, it is
%             indecisive, assign NaN 
%             FullClust(time).Predicted = NaN;
%             FullClust(time).CauseNaN.Edge = 1;
%         end
%     3.) If the cluster ID of the next time is NaN but the previous
%     time exists,...
%     elseif (isnan(FullClust(time+1).FltCluster))&&(~isnan(FullClust(time-1).FltCluster))
%         Look at the odd time segment after the next time
%         if FullClust(time-1).FltCluster == FullClust(time+2).FltCluster  % If the previous time = next 2 time
%            FullClust(time).Predicted = FullClust(time-1).FltCluster; % time-1 cluster is assigned in previous iteration
%            FullClust(time).CauseNaN.Edge = 0;
%         else
%             FullClust(time).Predicted = NaN;
%             FullClust(time).CauseNaN.Edge = 1; % Time segment during transition is indecisive
%         end
% 
%     else 
%          FullClust(time).Predicted = NaN;
%          FullClust(time).CauseNaN.Edge = 1; % Time segment during transition is indecisive
%     end
% end
% 
% 
% Plot time course of predicted clustering
% figure(4)
% stem([FullClust.Timeseg],[FullClust.Predicted],'Marker','.')
% axis([0 790 0 nclust+1])
% title('Clustering time course (predicted)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% 
% 4. Using moving window mode to smooth out the clustering
% WindowSize = 5;
% WindowEdge = floor(WindowSize/2);
% for time = 1:length(clustID)
%     Condition 1: Edges time segment where the candidate time does not
%     have 2 prior and 2 following time segments
%     if (time <= 2)&&(isnan(FullClust(time).Predicted)) % The starting time segments
%        window = [FullClust(1:WindowSize).Predicted];
%     elseif (time>=length(clustID)-1)&&(isnan(FullClust(time).Predicted))
%        window = [FullClust(length(clustID)-WindowSize:length(clustID)).Predicted];
%     elseif (time>=3) && (time<=length(clustID)-2)&&(isnan(FullClust(time).Predicted)) % Only trying to fill in NaN
%     Condition 2: The candidate is the middle one in the window, start at
%     3, ends at length(clustID)-2)
%        window = [FullClust(time-WindowEdge:time+WindowEdge).Predicted];
%     else
%        window = FullClust(time).Predicted;
%     end
%     FullClust(time).Smooth = mode(window);
% end
% 
% 
% Plot time course of predicted clustering
% figure(5)
% stem([FullClust.Timeseg],[FullClust.Smooth],'Marker','.','Color','b')
% axis([0 790 0 nclust+1])
% title('Clustering time course (predicted and smoothing)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% 
% Plot actual and predicted of cluster of one part
% figure(6)
% stem([FullClust(65:110).Timeseg],[FullClust(65:110).OrgCluster],'Marker','.','Color','b')
% hold on
% stem([FullClust(65:110).Timeseg],[FullClust(65:110).Smooth],'Marker','.','Color','r')
% axis([65 110 0 nclust+1])
% title('Clustering time course - comparison betweeb actual and predicted')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% legend('Actual cluster after artifact removal','Predicted cluster')
% hold off 
% 
% 6. Compare the predicted clustering with the actual clustering
% Check the accuracy/quality of expected and actual clustering
% 
% Inspector: Compare between predicted and actual cluster
% Inspector.correct = 0;    % If the expected and actual cluster match 
% Inspector.incorrect = 0;  % If the clusters are different and neither are NaN
% Inspector.nantime = 0;    % If either expected or actual is NaN.
% 
% for check = 1:length(clustID)
%     if (isnan(FullClust(check).Predicted)) || (isnan(FullClust(check).FltCluster))
%         Inspector.nantime = Inspector.nantime+1;    
%         FullClust(check).ExpectQuality = {'NaN'};   % Assign quality field
%     If actual cluster and expected cluster are the same,..
%     elseif FullClust(check).Predicted == FullClust(check).FltCluster 
%         Inspector.correct = Inspector.correct+1;    
%         FullClust(check).ExpectQuality = {'Consistent'};   % Assign quality field
%     If actual ad expected are different,'''
%     elseif FullClust(check).Predicted ~= FullClust(check).FltCluster
%         Inspector.incorrect = Inspector.incorrect+1;    
%         FullClust(check).ExpectQuality = {'Inconsistent'};   % Assign quality field
%     end
%  
% end
% 
% Calculate percentage of correct, incorrect and nantime for odd and even
% time
% Inspector.Total = length(clustID);
% 
% Percentage of each group
% Inspector.CorrectPer = 100*Inspector.correct/Inspector.Total;
% Inspector.IncorrectPer = 100*Inspector.incorrect/Inspector.Total;
% Inspector.NantimePer = 100*Inspector.nantime/Inspector.Total; %
% 
% 
% Plot the summary of consistency of clustering 
% Inspector.Summary = [Inspector.NantimePer, Inspector.IncorrectPer,Inspector.CorrectPer];
% Group={'NaN','Inconsistent','Consistent'};
% figure(7)
% bar(Inspector.Summary)
% set(gca,'XTickLabel',Group)
% Label value
% xpos = [1:3];
% ypos = [Inspector.Summary]+2;
% grouplabel = {num2str(Inspector.Summary(1),'%0.2f'),num2str(Inspector.Summary(2),'%0.2f'),num2str(Inspector.Summary(3),'%0.2f')};
% text(xpos,ypos,grouplabel,'HorizontalAlignment','center')
% title('Consistency of clustering with NaN data')
% Omit Nantime and calculate the consistency of non-nan cluster
% Inspector.TotalNN= Inspector.correct +Inspector.incorrect;
% Inspector.ConsistentPer = 100*Inspector.correct/Inspector.TotalNN;
% Inspector.InconsistentPer = 100*Inspector.incorrect/Inspector.TotalNN;
% Inspector.SummaryNN = [Inspector.InconsistentPer,Inspector.ConsistentPer];
% Consistency={'Inconsistent','Consistent'};
% figure(8)
% bar(Inspector.SummaryNN,'stacked')
% set(gca,'XTickLabel',Consistency)
% Label value
% xpos = [1:2];
% ypos = [Inspector.SummaryNN]+2;
% groupNNlabel = {num2str(Inspector.SummaryNN(1),'%0.2f'),num2str(Inspector.SummaryNN(2),'%0.2f')};
% text(xpos,ypos,groupNNlabel,'HorizontalAlignment','center')
% title('Consistency of clustering')
% 
% Alternative points (t_alt)
% o = 1; e = 1;
% for t_alt = 1:length(timedata)
%     if rem(t_alt,2)==0  % If time is even number
%        t_even(e) = t(t_alt);
%        e=e+1;
%     elseif rem(t_alt,2)~=0
%         t_odd(o) =t(t_alt);
%         o=o+1;
%     end
% end
% 
%  Extract odd time segment and even time segment for validation (?)
% Alternatively, create clustering vector of odd time and even time
% (PREVIOUSLY OBTAINED)
% 
% After obtaining the vectors of odd and even time, we can use odd time clustering to
% predict even and vice versa, but only for when there is no transition of
% stages (intermediate is undeterministic)
% t_predict = expected cluster number of t_even
% 
% Cluster ID at ODD time remains the same, analyse EVEN time
% [FullClust(t_odd).ExpectEven] = FullClust(t_odd).FltCluster; % Odd time and cluster ID remains the same
% 
% Cluster ID at EVEN time remains the same, analyse ODD time
% [FullClust(t_even).ExpectOdd] = FullClust(t_even).FltCluster; % Odd time and cluster ID remains the same
% 
%            *** Need a smarter fway for filtering and predicting***** 
% for time= 1:length(clustID)
%         If the ExpectOdd cluster ID is empty, i.e. ODD time (not assigned in previous step)
%     if isempty(FullClust(time).ExpectOdd)
%         Extra condition for first and last time segments
%         if (time==1)||(time==789)
%            FullClust(time).ExpectOdd = NaN;
%            FullClust(time).CauseNaN.Edge= 1;
%            continue;
%         end
%            
%         Conditions/ Rules for assigning cluster ID 
%         1.) If the two adjacent time segment have the same cluster ID, then
%         the middle one also have the same cluster
%         if FullClust(time-1).ExpectOdd == FullClust(time+1).ExpectOdd % If time is between time of stable stages (not during transition) 
%            FullClust(time).ExpectOdd = FullClust(time-1).ExpectOdd;
%            FullClust(time).CauseNaN.EvenTransition = 0;      % Transition between odd time
%         2.) If the cluster ID of previous time NaN but the cluster of
%         next time exists,...  
%         elseif (isnan(FullClust(time-1).ExpectOdd))&&(~isnan(FullClust(time+1).ExpectOdd)) %If the previous time segment is NaN, then look at the one before it
%             Look at the time segment before that to see if they are the
%             same.
%             if FullClust(time+1).ExpectOdd == FullClust(time-2).ExpectOdd
%                FullClust(time).ExpectOdd = FullClust(time-2).ExpectOdd; % time-2 cluster is assigned in previous iteration
%                FullClust(time).CauseNaN.EvenTransition = 0;
%             else
%                 If the time segment is at the transition, it is
%                 indecisive, assign NaN 
%                 FullClust(time).ExpectOdd = NaN;
%                 FullClust(time).CauseNaN.EvenTransition = 1;
%                 
%             end
%         3.) If the cluster ID of the next time is NaN but the previous
%         time exists,...
%         elseif (isnan(FullClust(time+1).ExpectOdd))&&(~isnan(FullClust(time-1).ExpectOdd))
%             Look at the odd time segment after the next time (cannot look
%             at time+2 as it has not been assigned yet)
%             if FullClust(time-1).ExpectOdd == FullClust(time+3).ExpectOdd  % If the previous time = next 2 time
%                FullClust(time).ExpectOdd = FullClust(time+3).ExpectOdd; % time-2 cluster is assigned in previous iteration
%                FullClust(time).CauseNaN.EvenTransition = 0;
%             else
%                 FullClust(time).ExpectOdd = NaN;
%                 FullClust(time).CauseNaN.EvenTransition = 1; % Time segment during transition is indecisive
%             end
% 
%         else 
%              FullClust(time).ExpectOdd = NaN;
%              FullClust(time).CauseNaN.EvenTransition = 1; % Time segment during transition is indecisive
%         end
%     end
%     
%     **************************************************************************************%
%     
%         If the ExpectEven cluster ID is empty, i.e. EVEN time (not assigned in previous step)
%     if isempty(FullClust(time).ExpectEven)
%         Conditions/ Rules for assigning cluster ID 
%         1.) If the two adjacent time segment have the same cluster ID, then
%         the middle one also have the same cluster
%         if FullClust(time-1).ExpectEven == FullClust(time+1).ExpectEven % If time is between time of stable stages (not during transition) 
%            FullClust(time).ExpectEven = FullClust(time-1).ExpectEven;
%            FullClust(time).CauseNaN.OddTransition = 0;      % Transition between odd time
%         2.) If the cluster ID of previous time NaN but the cluster of
%         next time exists,...  
%         elseif (isnan(FullClust(time-1).ExpectEven))&&(~isnan(FullClust(time+1).ExpectEven)) %If the previous time segment is NaN, then look at the one before it
%             Look at the time segment before that to see if they are the
%             same.
%             if FullClust(time+1).ExpectEven == FullClust(time-2).ExpectEven
%                FullClust(time).ExpectEven = FullClust(time-2).ExpectEven; % time-2 cluster is assigned in previous iteration
%                FullClust(time).CauseNaN.OddTransition = 0;
%             else
%                 If the time segment is at the transition, it is
%                 indecisive, assign NaN 
%                 FullClust(time).ExpectEven = NaN;
%                 FullClust(time).CauseNaN.OddTransition = 1;
%                 
%             end
%         3.) If the cluster ID of the next time is NaN but the previous
%         time exists,...
%         elseif (isnan(FullClust(time+1).ExpectEven))&&(~isnan(FullClust(time-1).ExpectEven))
%             Look at the odd time segment after the next time (cannot look
%             at time+2 as it has not been assigned yet)
%             if FullClust(time-1).ExpectEven == FullClust(time+3).ExpectEven  % If the previous time = next 2 time
%                FullClust(time).ExpectEven = FullClust(time+3).ExpectEven; % time-2 cluster is assigned in previous iteration
%                FullClust(time).CauseNaN.OddTransition = 0;
%             else
%                 FullClust(time).ExpectEven = NaN;
%                 FullClust(time).CauseNaN.OddTransition = 1; % Time segment during transition is indecisive
%             end
% 
%         else 
%              FullClust(time).ExpectEven = NaN;
%              FullClust(time).CauseNaN.OddTransition = 1; % Time segment during transition is indecisive
%         end
%     end
%     
% end
% 
% Plot time course of expect EVEN from ODD time clusters
% figure(6)
% stem([FullClust.Timeseg],[FullClust.ExpectEven],'Marker','.','Color','r')
% axis([0 790 0 nclust+1])
% title('Clustering time course (predicted EVEN, actual ODD)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% 
% Plot time course of expect ODD from EVEN time clusters
% figure(7)
% stem([FullClust.Timeseg],[FullClust.ExpectOdd],'Marker','.','Color','b')
% axis([0 790 0 nclust+1])
% title('Clustering time course (predicted ODD, actual EVEN)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% 
% 6. Compare the expected clustering with the actual clustering for both odd and even time 
% Check the accuracy/quality of expected and actual clustering
% 
% InspectEven: Compare between expected even time and even time clustering
% InspectEven.correct = 0;    % If the expected and actual cluster match 
% InspectEven.incorrect = 0;  % If the clusters are different and neither are NaN
% InspectEven.nantime = 0;    % If either expected or actual is NaN.
% 
% InspectOdd: Compare between expected odd time and odd time clustering
% InspectOdd.correct = 0;    % If the expected and actual cluster match 
% InspectOdd.incorrect = 0;  % If the clusters are different and neither are NaN
% InspectOdd.nantime = 0;    % If either expected or actual is NaN.
% 
% for check = 1:length(clustID)
%     2 cases: odd/even
%     switcher = rem(check,2);
%     switch switcher
%         case 0 % Checking EVEN time
%             If either expected and actual is NaN, the quality is
%             indecisive
%             if (isnan(FullClust(check).ExpectEven)) || (isnan(FullClust(check).FltCluster))
%                 InspectEven.nantime = InspectEven.nantime+1;    
%                 FullClust(check).ExpectEvenQuality = {'NaN'};   % Assign quality field
%             If actual cluster and expected cluster are the same,..
%             elseif FullClust(check).ExpectEven == FullClust(check).FltCluster 
%                 InspectEven.correct = InspectEven.correct+1;    
%                 FullClust(check).ExpectEvenQuality = {'Consistent'};   % Assign quality field
%             If actual ad expected are different,'''
%             elseif FullClust(check).ExpectEven ~= FullClust(check).FltCluster
%                 InspectEven.incorrect = InspectEven.incorrect+1;    
%                 FullClust(check).ExpectEvenQuality = {'Inconsistent'};   % Assign quality field
%             end
%         case 1 % Checking ODD time
%              If either expected and actual is NaN, the quality is
%             indecisive
%             if (isnan(FullClust(check).ExpectOdd)) || (isnan(FullClust(check).FltCluster))
%                 InspectOdd.nantime = InspectOdd.nantime+1;    
%                 FullClust(check).ExpectOddQuality = {'NaN'};   % Assign quality field
%             If actual cluster and expected cluster are the same,..
%             elseif FullClust(check).ExpectOdd == FullClust(check).FltCluster 
%                 InspectOdd.correct = InspectOdd.correct+1;    
%                 FullClust(check).ExpectOddQuality = {'Consistent'};   % Assign quality field
%             If actual ad expected are different,'''
%             elseif FullClust(check).ExpectOdd ~= FullClust(check).FltCluster
%                 InspectOdd.incorrect = InspectOdd.incorrect+1;    
%                 FullClust(check).ExpectOddQuality = {'Inconsistent'};   % Assign quality field
%             end
%     end
% end
% 
% Calculate percentage of correct, incorrect and nantime for odd and even
% time
% InspectEven.Total = length(t_even);
% InspectOdd.Total = length(t_odd);
% 
% Percentage of even time
% InspectEven.CorrectPer = InspectEven.correct/InspectEven.Total;
% InspectEven.IncorrectPer = InspectEven.incorrect/InspectEven.Total;
% InspectEven.NantimePer=InspectEven.nantime/InspectEven.Total; %
% 
% 
% = InspectEven.nantime/InspectEven.Total;
% 
% Percentage of odd time
% InspectOdd.CorrectPer = InspectOdd.correct/InspectOdd.Total;
% InspectOdd.IncorrectPer = InspectOdd.incorrect/InspectOdd.Total;
% InspectOdd.NantimePer = InspectOdd.nantime/InspectOdd.Total;
% 
% They should add up to 1; 
% tempOdd = InspectOdd.CorrectPer+InspectOdd.IncorrectPer+InspectOdd.NantimePer
% tempEven = InspectEven.CorrectPer+InspectEven.IncorrectPer+InspectEven.NantimePer
% 
% 
% Plot the correctness of each time segment
% Time course clustering with color-coded quality
% 
% Go through each time point and check the quality + cause of NaN and group
% them accordingly (to be plotted in different colour)
% Inspect EVEN (actual odd)
% n_OddAct=0;
% n_OddCon=0;
% n_OddInc=0;
% n_OddNaN=0;
% 
% n_EveAct=0;
% n_EveCon=0;
% n_EveInc=0;
% n_EveNaN=0;
% 
% Either group them now or assign group value before => another version
% 
% ***CAN BE COMBINED WITH FOR-LOOP ABOVE BUT SEPARATED NOT FOR CLEARER PURPOSE***
% for check = 1:length(clustID)
%     2 cases: odd/even 
%     switcher = rem(check,2);
%     switch switcher
%         case 0 % Checking EVEN time
%             If the check time is even, it is in ACTUAL group for
%             InspectOdd
%             n_OddAct=n_OddAct+1;
%             InspectOdd.ActualEven(n_OddAct) = check;
%             ============================= % 
%             If the check time is even, it is grouped in InspectEven according to quality 
%             if strcmp(FullClust(check).ExpectEvenQuality,'Consistent') % If the quality is 'Consistent'
%                 n_EveCon=n_EveCon+1;
%                 InspectEven.Consistent(n_EveCon)=check;
%             elseif strcmp(FullClust(check).ExpectEvenQuality,'Inconsistent') % If the quality is 'Consistent'
%                 n_EveInc=n_EveInc+1;
%                 InspectEven.Inconsistent(n_EveInc)=check;
%             elseif strcmp(FullClust(check).ExpectEvenQuality,'NaN') % If the quality is 'Consistent'
%                 n_EveNaN=n_EveNaN+1;
%                 InspectEven.NaNGroup(n_EveNaN)=check;
%             end
%         case 1 % Checking ODD time
%             If the check time is odd, it is in ACTUAL group for
%             InspectEven
%             n_EveAct=n_EveAct+1;
%             InspectEven.ActualOdd(n_EveAct) = check;
%              ============================= % 
%             If the check time is odd, it is grouped in InspectOdd according to quality 
%             if strcmp(FullClust(check).ExpectOddQuality,'Consistent') % If the quality is 'Consistent'
%                 n_OddCon=n_OddCon+1;
%                 InspectOdd.Consistent(n_OddCon)=check;
%             elseif strcmp(FullClust(check).ExpectOddQuality,'Inconsistent') % If the quality is 'Consistent'
%                 n_OddInc=n_OddInc+1;
%                 InspectOdd.Inconsistent(n_OddInc)=check;
%             elseif strcmp(FullClust(check).ExpectOddQuality,'NaN') % If the quality is 'Consistent'
%                 n_OddNaN=n_OddNaN+1;
%                 InspectOdd.NaNGroup(n_OddNaN)=check;
%                 Can expand this further for class of nan group
%             end
%     end
% end
% 
%  Once the time segment are grouped according to the quality, they are
% plotted with different colour
% Expect EVEN, Actual ODD
% figure(8)
% Actual value 
% stem([InspectEven.ActualOdd],[FullClust(InspectEven.ActualOdd).OrgCluster],'Color','k','Marker','.')
% axis([0 790 0 nclust+1])
% title('Clustering time course (with quality of even time)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% hold on
% Consisten with adjacent time
% stem([InspectEven.Consistent],[FullClust(InspectEven.Consistent).OrgCluster],'Color','r','Marker','.')
% Inconsistent with adjacent time
% stem([InspectEven.Inconsistent],[FullClust(InspectEven.Inconsistent).OrgCluster],'Color','b','Marker','x')
% NaN due to artifact filtering, edge time or during transition
% stem([InspectEven.NaNGroup],[FullClust(InspectEven.NaNGroup).OrgCluster],'Color','g','Marker','.')
% legend('Actual ODD time cluster','Consistent','Inconsistent','NaN')%,'Location','NorthEastOutside')
% hold off
% 
%  ===============================================
% Expect ODD, Actual EVEN
% figure(9)
% Actual value 
% stem([InspectOdd.ActualEven],[FullClust(InspectOdd.ActualEven).OrgCluster],'Color','k','Marker','.')
% axis([0 790 0 nclust+1])
% title('Clustering time course (with quality of odd time)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% hold on
% Consisten with adjacent time
% stem([InspectOdd.Consistent],[FullClust(InspectOdd.Consistent).OrgCluster],'Color','r','Marker','.')
% Inconsistent with adjacent time
% stem([InspectOdd.Inconsistent],[FullClust(InspectOdd.Inconsistent).OrgCluster],'Color','b','Marker','x')
% NaN due to artifact filtering, edge time or during transition
% stem([InspectOdd.NaNGroup],[FullClust(InspectOdd.NaNGroup).OrgCluster],'Color','g','Marker','.')
% legend('Actual EVEN time cluster','Consistent','Inconsistent','NaN')%,'Location','NorthEastOutside')
% hold off
% 
% Combine the quality plot of ODD and EVEN time
% Similar to above but omit actual even/odd
% figure(10)
% axis([0 790 0 nclust+1])
% title('Clustering time course (with quality of ALL time)')
% xlabel('Epoch (30-second segment)')
% ylabel('Cluster #')
% hold on
% Consisten with adjacent time
% stem([InspectOdd.Consistent],[FullClust(InspectOdd.Consistent).OrgCluster],'Color','r','Marker','.')
% Inconsistent with adjacent time
% stem([InspectOdd.Inconsistent],[FullClust(InspectOdd.Inconsistent).OrgCluster],'Color','b','Marker','x')
% NaN due to artifact filtering, edge time or during transition
% stem([InspectOdd.NaNGroup],[FullClust(InspectOdd.NaNGroup).OrgCluster],'Color','g','Marker','.')
% legend('Consistent','Inconsistent','NaN')%,'Location','NorthEastOutside')
% Consisten with adjacent time
% stem([InspectEven.Consistent],[FullClust(InspectEven.Consistent).OrgCluster],'Color','r','Marker','.')
% Inconsistent with adjacent time
% stem([InspectEven.Inconsistent],[FullClust(InspectEven.Inconsistent).OrgCluster],'Color','b','Marker','x')
% NaN due to artifact filtering, edge time or during transition
% stem([InspectEven.NaNGroup],[FullClust(InspectEven.NaNGroup).OrgCluster],'Color','g','Marker','.')
% hold off