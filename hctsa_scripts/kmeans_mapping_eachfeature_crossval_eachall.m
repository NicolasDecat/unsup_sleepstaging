% Changes from original script
% - Changed configuration settings
% - epochSelect (instead of epochSelectFunction)
% - for ch=1:NUM_CHANNELS_TO_RUN, instead of ch=1:number_of_channels_used
% - clear max (because was used as a variable beforehand)

%% Function: Count the number of epochs in each stages and recore the epochIDs

function [testMat Nf testTS PERC_STAGE] = cross_validation_selectivefeatures_crossval_eachall(experiment, hctsa_ops, cm_save_dir, number_of_channels_used, epochSelectFunction, selective_feature,sub,v,col)

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
            trainMat = hctsa_ops(trainTS',:);
            testMat{Nf} = hctsa_ops(testTS',:);
        else
            % increment=(size(hctsa_ops,1)/3)*(ch-1);
            increment=round((size(hctsa_ops,1)/7))*(ch-1);
            trainMat = [trainMat hctsa_ops((trainTS+increment)',:)];
            testMat(:,Nf) = [testMat hctsa_ops((testTS+increment)',:)];
        end
    end
    
    
  
end % End Nf-th randomisation


%% Compute type1AUC

   run('type1auc_eachfeature_crossval_eachall.m')


end