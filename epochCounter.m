%% Function: Count the number of epochs in each stages and recore the epochIDs
% Input: - whichData // which dataset to be used ([1,5,7,13,14] for now)
%        - LabeledStage // labelled sleep stage from annotation
% Output : stgID // randomised order of epoch IDs

function stgID = epochCounter(whichData,LabeledStage)
%% Check if data set exists
% Use ccshs1800001 as default (avoid error)
validData = [1,5,7,13,14]; % Valid datasets
if ~ismember(whichData,validData)
    % whichData = 1; % Default case 
    % *** can change to display error and keep asking for valid input using
    % Add while loop?
    error('Dataset is invalid...input: 1,5,7,13,14\n')
end
%% Remove initial W stage from randomisation and sampling
% Marking the end of W stage and Sleep stages (N1-N3) (Manually)
endW = [334,380,391,375,174];
endS = [1374,1442,1442,1531,1492];
endID = find(whichData==validData); % Determine which ending/beginning to use

selectID = endW(endID)+1:endS(endID)-1;

% Selected ID and their labels
selectLabel = LabeledStage(selectID);

% Record in struct
stgID.selectID = selectID;
stgID.selectLabel = selectLabel;

%% Counting number of stages
% Proportion of each sleep stage (0 - wake, 1-4 NREM, 5 - REM)
stgNum = length(unique(LabeledStage));
stgLab = {'W','N1','N2','N3','R'};

% Record number of epochs in each stage and the ID of the epoch 
n = zeros(stgNum,1);
for i = 1:length(selectLabel)
    for j = 1:stgNum
        if selectLabel(i)<5
            if selectLabel(i)==(j-1) % Stage 0-3 (W,N1,N2,N3)
                n(j)=n(j)+1;    % Increment number of epoch at this stage
                stgID.allID.(stgLab{j})(n(j))=i;    % Record epoch IDs
            end
        else
            if selectLabel(i)== j   % Stage 5 (R)
                n(j)=n(j)+1;    % Increment number of epochs
                stgID.allID.(stgLab{j})(n(j))=i;
            end
        end
    end
end
stgID.nStg = n;

%% Remove class with less than cut-off
% cutoff = 0.02*length(LabeledStage);   % Exclude from the analysis if the stage occur less than 2% of total sleep period
% n=0;
% for i = 1:length(stgID.nStg)
%     if stgID.nStg(i)<cutoff
% %         n=n+1;
% %         stgID.useStg(n)=i;
% %         stgID.usePro(n)=stgID.stgPro(i);
% %         stgID.useStgName(n) = stgLab(i);
%     end
% end

%% Minimum samples
stgID.Nmin = min(stgID.nStg);

%% \Randomise order of each stage/
% ********** Uncomment this section if crossvalidation requires same epochs for
% each run (Only re-sample training and test datasets)
% If different sets of epochs are used in each iteration of the loop use
% code in epochSelect.m *******************

% % Randomised epochID of each class is in stgID.actualID
% % these epochIDs will be used for training and testing
% stgID.useID =[];
% 
% for m=1:stgNum
%     randID = randperm(stgID.nStg(m),stgID.Nmin); % Random permutation - index of one stage
%     allID = stgID.allID.(stgLab{m}); % All epoch id of that stage
%     useID = allID(randID);  % Take randomised epoch id
%     stgID.useID.(stgLab{m}) = useID;    % Index in selectLabel
%     stgID.actualID.(stgLab{m}) = selectID(useID);
% end
end