% Changes from original script
% - Changed configuration settings
% - epochSelect (instead of epochSelectFunction)
% - for ch=1:NUM_CHANNELS_TO_RUN, instead of ch=1:number_of_channels_used
% - clear max (because was used as a variable beforehand)

%% Function: Count the number of epochs in each stages and recore the epochIDs

function [testMat Nf testTS trainTS trainMat label stats PERC_STAGE] = cross_validation_selectivefeatures_crossval_eachall_SVM(experiment, hctsa_ops, cm_save_dir, number_of_channels_used, epochSelectFunction, selective_feature,sub,v,col,D)

%% Cross-validation code

configuration_settings;


%% Obtain epochID using epochCounter() function
whichData = WHICH_DATA;
annotation = load(ANSWER_FILE);
label = annotation.sleepstage;
stgID = epochCounter(whichData,label);
stgLab = {'W','N1','N2','N3','R'};

% Training
trainingProportion = TRAINING_PERCENTAGE;
nIterations = 10;

%% Multiple iteration of randomisation and cross-validation
% Initialise result struct
block(nIterations) = struct();
stats = struct();

for Nf = 1:nIterations
    
    [block(Nf).trainTS,block(Nf).testTS]=epochSelect(stgID,trainingProportion);
    
    % The following turn the nxm matrix to 1x(n*m) matrix
    trainTS = block(Nf).trainTS.';
    trainTS = trainTS(:).';
    testTS = block(Nf).testTS.';
    testTS = testTS(:).';
  
    
    %% Select data of wanted time ID  
  
    NUM_CHANNELS_TO_RUN = v;
    
    for ch=1:NUM_CHANNELS_TO_RUN

        if ch==1
            trainMat{Nf} = hctsa_ops(trainTS',:);
            testMat{Nf} = hctsa_ops(testTS',:);
        else
            % increment=(size(hctsa_ops,1)/3)*(ch-1);
            increment=round((size(hctsa_ops,1)/7))*(ch-1);
            trainMat = [trainMat hctsa_ops((trainTS+increment)',:)];
            testMat(:,Nf) = [testMat hctsa_ops((testTS+increment)',:)];
        end
    end
    
    stats.scoredTrain(Nf,:) = label(trainTS)';
    stats.scoredTest(Nf,:) = label(testTS)';
    
    %% Reconstruct testMat and remove SV features
    onlyWB = true;
    
    if onlyWB == true
            load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/allfeat_removed')
            load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ALL_removed_feat(2146)')

            TRY_test = testMat{1,Nf};
            TRY_train = trainMat{1,Nf};

            InsertCol = zeros(size(TRY_test,1),1);
            InsertCol_train = zeros(size(TRY_train,1),1);


            for x = 1:length(removed_feat_idx{1,D})

                Part1 = TRY_test(:,1:removed_feat_idx{1,D}(x)-1);
                Part2 = [InsertCol TRY_test(:,removed_feat_idx{1,D}(x):end)];
                TRY_test = ([Part1 Part2]);
                removed_feat_idx{1,D} = removed_feat_idx{1,D};

                Part1_train = TRY_train(:,1:removed_feat_idx{1,D}(x)-1);
                Part2_train = [InsertCol_train TRY_train(:,removed_feat_idx{1,D}(x):end)];
                TRY_train = ([Part1_train Part2_train]);
                removed_feat_idx{1,D} = removed_feat_idx{1,D};

            end

            % Remove both the specific and commonly removed features
            TRY_test(:,spec_and_common_feat) = []; 
            TRY_train(:,spec_and_common_feat) = []; 

            testMat(1,Nf) = {TRY_test};
            trainMat(1,Nf) = {TRY_train};
    end
  
end % End Nf-th randomisation


%% Compute type1AUC

   run('type1auc_eachfeature_crossval_eachall_SVM.m')


end