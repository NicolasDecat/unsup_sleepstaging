%%%%%
%%%%%

%% Type 1 AUC hctsa for each feature: one-vs-one paired classification

if isfile('/Users/nico/Documents/HCTSA/Analysis/AUC/Per_correct_mean.mat') == 1
    load('/Users/nico/Documents/HCTSA/Analysis/AUC/Per_correct_mean.mat')
end

IsStage1 = [];
IsStage2 = [];

% Find how many epochs for each sleep stage
stgL = size(testMat,1);
stgL = stgL/5;          

% Get epochs indices for each stage
wake = (1:stgL);
N1 = ((stgL+1):(stgL*2));
N2 = ((stgL*2+1):(stgL*3));
N3 = ((stgL*3+1):(stgL*4));
rem = ((stgL*4+1):(stgL*5));

% All binary classifiers
allpairs = [{[wake;N1]} {[wake;N2]} {[wake;N3]} {[wake;rem]} {[N1;N2]} {[N1;N3]} {[N1;rem]} {[N2;N3]} {[N2;rem]} {[N3;rem]}];

for Nit = 1:size(testMat,2)
    
    testMatr = testMat(:,Nit);  % One iteration of testMat at a time

    for C = 1:10         % For each binary classifier
        
        % Choose one classifier
        classifier = allpairs{C};     

        for test = 1:stgL           % each epoch will be taken out and used as testing epoch
        
            % Get epoch indices for both stages of the classifier
            stage1 = classifier(1,:);    
            stage2 = classifier(2,:);   
            
            % Testing data (1 epoch for each stage is taken out: "leave-1-out" strategy for cross validation)
            test_stage1 = stage1(test);   % The only piece of information that we change for each 'test' iteration
            test_stage2 = stage2(test);
 
            % Training data (N-1 epochs for each stage)
            stage1 = setdiff(stage1,test_stage1);   % Leave one epoch out ("test" epoch)
            stage2 = setdiff(stage2,test_stage2);  
            
            % Get median of hctsa response for the two stages
            hctsa_resp1 = testMatr(stage1,1);
            hctsa_resp2 = testMatr(stage2,1);
         
            hctsa_stage1 = median(hctsa_resp1);   
            hctsa_stage2 = median(hctsa_resp2);

            [~,index_max] = max([hctsa_stage1 hctsa_stage2]);   % Will be useful later to set direction
            
            % Get the threshold: value between the 2 medians
            thresh = (hctsa_stage1 + hctsa_stage2)/2;   

            % Assign testing epoch of Stage 1 to Stage 1 or 2
            if testMatr(test_stage1) > thresh && index_max == 1     % If hctsa response of testing epoch is above threshold and if mean hctsa response of training epochs from Stage 1 is higher than Stage 2                          % If hctsa response of training epoch is the highest value, then being above threshold means bleonging to stage1
                   IsStage1 = [IsStage1 test_stage1];
            elseif testMatr(test_stage1) > thresh && index_max == 2
                   IsStage2 = [IsStage2 test_stage1];
            elseif testMatr(test_stage1) < thresh && index_max == 1   
                   IsStage2 = [IsStage2 test_stage1];
            elseif testMatr(test_stage1) < thresh && index_max == 2   
                   IsStage1 = [IsStage1 test_stage1];
            end

            % Assign testing epoch of Stage 2 to Stage 1 or 2
            if testMatr(test_stage2) > thresh && index_max == 1     % If hctsa response of testing epoch is above threshold and if mean hctsa response of trainning epochs from Stage 1 is higher than Stage 2                          % If hctsa response of training epoch is the highest value, then being above threshold means bleonging to stage1
                   IsStage1 = [IsStage1 test_stage2];
            elseif testMatr(test_stage2) > thresh && index_max == 2
                   IsStage2 = [IsStage2 test_stage2];
            elseif testMatr(test_stage2) < thresh && index_max == 1   
                   IsStage2 = [IsStage2 test_stage2];
            elseif testMatr(test_stage2) < thresh && index_max == 2   
                   IsStage1 = [IsStage1 test_stage2];
            end

            %%% Store assigned stages (to double check later if needed)
%             hctsa_stage1_stored(C) = hctsa_stage1;
%             hctsa_stage2_stored(C) = hctsa_stage2;
%             thresh_stored(C) = thresh;
            
%             IsStage1_stored{C} = IsStage1;
%             IsStage2_stored{C} = IsStage2;
            
            
        end

         % 1 binary classifier is done: calculate % correct
         Per_correct_stage1(Nit,C) = (length(find(ismember(IsStage1,classifier(1,:))))/stgL)*100;   % Percentage of epochs from stage 1 labeled right
         Per_correct_stage2(Nit,C) = (length(find(ismember(IsStage2,classifier(2,:))))/stgL)*100;   % Percentage of epochs from stage 1 labeled right 

         IsStage1 = [];
         IsStage2 = [];

    end
    
    % Mean across all iterations
    Per_correct_stage1_mean = mean(Per_correct_stage1);
    Per_correct_stage2_mean = mean(Per_correct_stage1);
    
    
end  

% Mean across iterations across stage 1 and 2 (not sure if correct), % for 1 feature
Per_correct_mean(:,FF) = mean([Per_correct_stage1_mean;Per_correct_stage2_mean])';
 
% % Save
fpath = '/Users/nico/Documents/HCTSA/Analysis/AUC';
save(fullfile(fpath,'Per_correct_mean.mat'),'Per_correct_mean')  % Save all columns in AUC folder


%% Plot distribution curves

% figure; h = histfit(hctsa_resp1); hold on; j = histfit(hctsa_resp2);
% 
% hold on
% median_stage_1 = hctsa_stage1;
% xline(median_stage_1,'b--','LineWidth',2)
% median_stage_2 = hctsa_stage2;
% xline(median_stage_2,'g--','LineWidth',2)
% threshold = thresh;
% xline(threshold,'k--','LineWidth',2)
% 
% % Color stage 1
% h(1).FaceColor = [.9 .9 1]; % histogram (light blue)
% h(2).Color = [0 0 1];  % curve (blue)
% % Color stage 2
% j(1).FaceColor = [.9 1 .9]; % histogram (light green)
% j(2).Color = [0 1 0];  % curve (green)
% 
% hold off
% 
% disp('ok')
