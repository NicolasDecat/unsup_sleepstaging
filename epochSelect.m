%% Function: Divide epochs into 70% training and 30% test sets
% After running epochSampX() outsie Nf loop (epochSampX has to be called
% only once for each dataset)
% Input: stgID (output of epochSampX), proportion of training data [0-1]
% Output: trainID, testID (actual ID of epochs to be used when selecting
% data for cross-validation anaylsis)

% function [trainTS, testTS] = epochSelect(stgID,pTraining)
function [trainTS, testTS] = epochSelect(stgID,trainingProportion)

% Stage labels
stgLab = {'W','N1','N2','N3','R'};
stgNum = length(unique(stgID.selectLabel));

% trainPP = round(pTraining*stgID.Nmin); 
trainPP = round(trainingProportion*stgID.Nmin); 

%% Randomly sample from each class
for m=1:stgNum
    %% \Randomise order of each stage/
    % ******* Uncomment this section if same set of epochs are to be used
    % in every iteration of cross validation. **************** 
    randID = randperm(stgID.nStg(m),stgID.Nmin);    % Random permutation - index of one stage: We randomly generate 37 (Nmin) numbers from 1 to 242 (nb of W epochs): indexes saved as randID 
    allID = stgID.allID.(stgLab{m});                % All possible epoch IDs of the stage
    useID = allID(randID);                          % Take randomised epoch id: the 37 random numbers (randID) pick inside all possible W epochs (allID)
    stgID.useID.(stgLab{m}) = useID;                % the randomly selected W epochs saved as useID (useID = random IDs from W total epochs)
    stgID.actualID.(stgLab{m}) = stgID.selectID(useID);   % Actual epochIDs   (actual IDs = equivalent W IDs when picking amont total recording epochs.)
    
    %% Sample 70% for training and 30% for test
    sampID = stgID.useID.(stgLab{m}); % Redundant but in the case of switching between two versions
    randtrain = randperm(stgID.Nmin);
    trainID = sampID(randtrain(1:trainPP)); % Select actual ID as 70% training: take 26/37 from useID (the randomly selected W epochs)
    testID = sampID(randtrain((trainPP+1):end));    % Select actual epoch IDs as test
    stgID.trainID.(stgLab{m})= trainID;
    stgID.testID.(stgLab{m})= testID;
    % Combine training ID to perform classification
    if m==1 % First iteration
        trainTS = stgID.selectID(trainID);  
        testTS = stgID.selectID(testID);
    else
        trainTS = [trainTS;stgID.selectID(trainID)]; 
        testTS = [testTS;stgID.selectID(testID)];  % join all random testing epochs from all 5 stages
    end
end
end