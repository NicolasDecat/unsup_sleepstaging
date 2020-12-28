
%% Type 1 AUC hctsa for each feature: one-vs-one paired classification


if isfile('/Users/nico/Documents/HCTSA/Analysis/AUC/AUC_per_feature.mat') == 1
    load('/Users/nico/Documents/HCTSA/Analysis/AUC/AUC_per_feature.mat')
else
    AUC_per_feature = zeros(10,10,1);  % 10 classifiers x 10 features x 1 channel conditions
end

AUC = zeros(10,1);


for Nit = 1:size(scoredTest,1)   % For each iteration Nf

      
    ORIGlabels = scoredTest(Nit,:);  
    CLUSTERdecision = predictTest(Nit,:);

    % Number of epochs for each sleep stage
    stgL = size(testMat,1);
    stgL = stgL/5;              
    
    % Indexes for each sleep stage
    wake = (1:stgL);
    N1 = ((stgL+1):(stgL*2));
    N2 = ((stgL*2+1):(stgL*3));
    N3 = ((stgL*3+1):(stgL*4));
    rem = ((stgL*4+1):(stgL*5));

    % 0 vs 1
    origlabels = ORIGlabels(1,[wake,N1]);
    clusterdecision = CLUSTERdecision(1,[wake,N1]);
    posclass = 1;    
    [X,Y,T,AUC(1,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 0 vs 2
    origlabels = ORIGlabels(1,[wake,N2]);
    clusterdecision = CLUSTERdecision(1,[wake,N2]);
    posclass = 2;    
    [X,Y,T,AUC(2,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 0 vs 3
    origlabels = ORIGlabels(1,[wake,N3]);
    clusterdecision = CLUSTERdecision(1,[wake,N3]);
    posclass = 3;    
    [X,Y,T,AUC(3,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 0 vs rem
    origlabels = ORIGlabels(1,[wake,rem]);
    clusterdecision = CLUSTERdecision(1,[wake,rem]);
    posclass = 5;    
    [X,Y,T,AUC(4,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 1 vs 2
    origlabels = ORIGlabels(1,[N1,N2]);
    clusterdecision = CLUSTERdecision(1,[N1,N2]);
    posclass = 2;    
    [X,Y,T,AUC(5,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 1 vs 3
    origlabels = ORIGlabels(1,[N1,N3]);
    clusterdecision = CLUSTERdecision(1,[N1,N3]);
    posclass = 3;    
    [X,Y,T,AUC(6,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 1 vs rem
    origlabels = ORIGlabels(1,[N1,rem]);
    clusterdecision = CLUSTERdecision(1,[N1,rem]);
    posclass = 5;    
    [X,Y,T,AUC(7,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 2 vs 3
    origlabels = ORIGlabels(1,[N2,N3]);
    clusterdecision = CLUSTERdecision(1,[N2,N3]);
    posclass = 3;    
    [X,Y,T,AUC(8,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 2 vs rem
    origlabels = ORIGlabels(1,[N2,rem]);
    clusterdecision = CLUSTERdecision(1,[N2,rem]);
    posclass = 5;    
    [X,Y,T,AUC(9,Nit)] = perfcurve(origlabels,clusterdecision,posclass);

    % 3 vs rem
    origlabels = ORIGlabels(1,[N3,rem]);
    clusterdecision = CLUSTERdecision(1,[N3,rem]);
    posclass = 5;    
    [X,Y,T,AUC(10,Nit)] = perfcurve(origlabels,clusterdecision,posclass);
    
end


% Get mean AUC for each classifier
for rowAUC = 1:10                        % For all rows of AUC (all 10 classifiers)
    Mean_AUC(rowAUC,FF) = mean(AUC(rowAUC,:));
end

% Save

AUC_per_feature(1:10,FF,v) = Mean_AUC(:,FF); % Rows = binary classifiers in order listed above (row 1 = 0vs1, row 2 = 0vs2, etc)
fpath = '/Users/nico/Documents/HCTSA/Analysis/AUC';
save(fullfile(fpath,'AUC_per_feature.mat'),'Mean_AUC','AUC_per_feature')  % Save all columns in AUC folder

% See A(:,:,v), and know that {AUC_per_feature} gives what you want




