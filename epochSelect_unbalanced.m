%% Function: Divide epochs into 70% training and 30% test sets
% After running epochSampX() outsie Nf loop (epochSampX has to be called
% only once for each dataset)
% Input: stgID (output of epochSampX), proportion of training data [0-1]
% Output: trainID, testID (actual ID of epochs to be used when selecting
% data for cross-validation anaylsis)

function [trainTS, testTS] = epochSelect(stgID,pTraining)
% Stage labels
stgLab = {'W','N1','N2','N3','R'};
stgNum = length(unique(stgID.selectLabel));

trainPP = round(pTraining*stgID.Nmin); 
%% Randomly sample from each class
for m=1:stgNum
    %% \Randomise order of each stage/
    % ******* Uncomment this section if same set of epochs are to be used
    % in every iteration of cross validation. **************** 
    randID = randperm(stgID.nStg(m),stgID.Nmin);    % Random permutation - index of one stage
    allID = stgID.allID.(stgLab{m});                % All epoch id of that stage
    useID = allID(randID);                          % Take randomised epoch id
    stgID.useID.(stgLab{m}) = useID;                % Index in selectLabel
    stgID.actualID.(stgLab{m}) = stgID.selectID(useID);   % Actual epochIDs
    %
    %% Sample 70% for training and 30% for test
    sampID = stgID.useID.(stgLab{m}); % Redundant but in the case of switching between two versions
    randtrain = randperm(stgID.Nmin);
    trainID = sampID(randtrain(1:trainPP)); % Select actual ID as 70% training
%     testID = sampID(randtrain((trainPP+1):end));    % Select actual epoch IDs as test
    stgID.trainID.(stgLab{m})= trainID;
    %stgID.testID.(stgLab{m})= testID;
    % Combine training ID to perform classification
    if m==1 % First iteration
        trainTS = stgID.selectID(trainID);
%         testTS = setdiff(stgID.selectID, trainID);
    else
        trainTS = [trainTS;stgID.selectID(trainID)]; 
    end
end
testTS = setdiff(stgID.selectID, trainTS);

end