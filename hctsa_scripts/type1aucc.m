
%% Type 1 AUC hctsa
 
if isfile('/Users/nico/Documents/HCTSA/Analysis/iterdata.mat') == 1
    load('/Users/nico/Documents/HCTSA/Analysis/iterdata.mat')
    Iteration = table2array(SummaryTable(:,1));
    NumChannels = table2array(SummaryTable(:,2));                  
    Dataset = table2array(SummaryTable(:,3));                                
    Sleep_stage = table2array(SummaryTable(:,4));                         
    Testing_accuracy = table2array(SummaryTable(:,5));
    AUC = table2array(SummaryTable(:,6)); 
else
    col = 1;
end

% 0: wake, 1: N1, 2: N2, 3: N3, 5: REM
idx = [0 1 2 3 5];


for Nit = 1:size(scoredTest,1)   % For each iteration Nf

    col_idx = 1;
    
    for i = idx(1:5)   % For all 5 stages
  
        origlabels = scoredTest(Nit,:);  
        clusterdecision = predictTest(Nit,:);  

        % New labels for original labels
        idxstg = find(origlabels==i); 
        idxrest = find(origlabels~=i);
        origlabels(idxstg) = 6;         % 6 = label given to selected stage
        origlabels(idxrest) = 7;        % 7 = label given to all other stages

        % New labels for cluster decision
        idxstg2 = find(clusterdecision==i); 
        idxrest2 = find(clusterdecision~=i);
        clusterdecision(idxstg2) = 6;
        clusterdecision(idxrest2) = 7;
        posclass = 7;

        % Compute AUC
        [X,Y,T,stgAUC(col,1)] = perfcurve(origlabels,clusterdecision,posclass);

        %Save
        Iteration(col,1) = Nit;                                             % Col 1: number of iterations
        Dataset(col,1) = convertCharsToStrings(sub);                        % Col 2: Dataset ID
        NumChannels(col,1) = NUM_CHANNELS_TO_RUN;                           % Col 3: Number of channel
        Sleep_stage(col,1) = idx(col_idx);                                  % Col 5: stage for ind AUC
        Testing_accuracy(col,1) = ((sum((scoredTest(Nit,:) == predictTest(Nit,:))'))/size(scoredTest, 2))';
        AUC(col,1) = stgAUC(col);                                           % Col 6: AUC                                  % Col 8: predict test (for conf matrix)
        col = col+1;
        col_idx = col_idx+1;

    end   
    
end

SummaryTable = table(Iteration,NumChannels,Dataset,Sleep_stage,Testing_accuracy,AUC);

fpath = '/Users/nico/Documents/HCTSA/Analysis';
save(fullfile(fpath,'iterdata.mat'),'SummaryTable','col')  % Save both variables in AUC folder

% Average of all AUC
% AUC(6,1) = zeros(1);
% AUC(6,1) = mean(AUC(1:5,1));
% AUC = AUC(6,1);


%%% Using R

% response = origlables(1,:)';  
% predictor = stats.predictTest(1,:)'; 


% response = origlabels';  
% predictor = clusterdecision'; 
% 
% Test = [zeros(size(response)) zeros(size(response))];
% Test(:,1) = response;
% Test(:,2) = predictor;
% Test = array2table(Test);
% Test.Properties.VariableNames{1} = 'response';
% Test.Properties.VariableNames{2} = 'predictor';
% 
% writetable(Test,'Test.xlsx')
% 

%% Manual one-vs-one 
% 
% stgL = size(testMat,1);
% stgL = stgL/5;
% 
% AUC = zeros(11,1);
% wake = (1:stgL);
% N1 = ((stgL+1):(stgL*2));
% N2 = ((stgL*2+1):(stgL*3));
% N3 = ((stgL*3+1):(stgL*4));
% rem = ((stgL*4+1):(stgL*5));
% 
%   
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 0 vs 1
% origlabels = origlabels(1,[wake,N1]);
% clusterdecision = clusterdecision(1,[wake,N1]);
% posclass = 1;    
% [X,Y,T,AUC(1,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 0 vs 2
% origlabels = origlabels(1,[wake,N2]);
% clusterdecision = clusterdecision(1,[wake,N2]);
% posclass = 2;    
% [X,Y,T,AUC(2,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 0 vs 3
% origlabels = origlabels(1,[wake,N3]);
% clusterdecision = clusterdecision(1,[wake,N3]);
% posclass = 3;    
% [X,Y,T,AUC(3,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 0 vs rem
% origlabels = origlabels(1,[wake,rem]);
% clusterdecision = clusterdecision(1,[wake,rem]);
% posclass = 5;    
% [X,Y,T,AUC(4,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 1 vs 2
% origlabels = origlabels(1,[N1,N2]);
% clusterdecision = clusterdecision(1,[N1,N2]);
% posclass = 2;    
% [X,Y,T,AUC(5,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 1 vs 3
% origlabels = origlabels(1,[N1,N3]);
% clusterdecision = clusterdecision(1,[N1,N3]);
% posclass = 3;    
% [X,Y,T,AUC(6,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 1 vs rem
% origlabels = origlabels(1,[N1,rem]);
% clusterdecision = clusterdecision(1,[N1,rem]);
% posclass = 5;    
% [X,Y,T,AUC(7,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 2 vs 3
% origlabels = origlabels(1,[N2,N3]);
% clusterdecision = clusterdecision(1,[N2,N3]);
% posclass = 3;    
% [X,Y,T,AUC(8,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 2 vs rem
% origlabels = origlabels(1,[N2,rem]);
% clusterdecision = clusterdecision(1,[N2,rem]);
% posclass = 5;    
% [X,Y,T,AUC(9,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% origlabels = scoredTest;  
% clusterdecision = predictTest;
% % 3 vs rem
% origlabels = origlabels(1,[N3,rem]);
% clusterdecision = clusterdecision(1,[N3,rem]);
% posclass = 5;    
% [X,Y,T,AUC(10,1)] = perfcurve(origlabels,clusterdecision,posclass);
% 
% AUC(11,1) = mean(AUC(1:10,1));
% AUC = AUC(11,1)
