function TransitionMat = TPM(DataState,nState)
% Transition probability matrix function
% Generate transition probability matrix of states/clusters
% 
% Input: state vector matrix of data, number of state
% Include code to check number of state 
% Output: Transition probability matrix
% Variables: DataState = state vetor matrix of the data, nState = number of
% state, TransitionMat = transition probability matrix 
% 
% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% Check if the number of states specified matched with the number of state
% in data state
if (max(DataState)>nState)
   fprintf('Number of states detected is more than specified. Adjusting number of state...\n')
elseif (max(DataState)<nState)
    fprintf('Some states are invalid. Adjusting number of state...\n')
end

% Finalise number of states
nState = max(DataState);
fprintf('This data set has %d states/clusters.\n',nState)

% Initialise TPM (n states = nxn matrix)
Tmat = zeros(nState,nState);

for m=1:length(DataState)-1         % Number of transition is one less than number of data points
    currentstate = DataState(m);    % Cluster number at current t (m=t)
    nextstate = DataState(m+1);     % Cluster number at the next t (m=t+1)
    
    Tmat(currentstate,nextstate) = Tmat(currentstate,nextstate)+1; % Increment the number of transition
end

% Normalise transition to get probability of transtion

% Number of data points in each of the current state/cluster
statesum = sum(Tmat,2);  % Horizontal sum (sum of each states)

% Obtain probability matrix 
normMatrix = repmat(statesum,1,nState);
TransitionMat = Tmat./normMatrix;

return