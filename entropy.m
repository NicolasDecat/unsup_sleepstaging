function H = entropy(x,option)
% H - entropy of data
% x is the transitional probability matrix
% ===========================
% Size of TPM - number of states
[current, next]=size(x);

% Check option
if isempty(option)
   option = 'default';
end

% Calculate entropy
if strcmp(option,'default')
    for i=1:current for j=1:next
            if x(i,j)~=0
                H_ind(i,j)=x(i,j)*log2(x(i,j));  % Value of each element
            else
                H_ind(i,j)=0;
            end
        end
    end
    % Entropy: 
    H = -sum(sum(H_ind));
    return
elseif strcmp(option,'conditional')
    for i=1:current   % For each current state, the probability of transitioning to next states is x(i,j)
        for j=1:next
            if x(i,j)~=0    % If probability is non-zero.
                H_ind(j)=x(i,j)*log2(x(i,j));
               % H_ind(i,j)=x(i,j)*log2(x(i,i)/x(i,j));  % Value of each element
            else
               H_ind(j)=0;
            end
        end
        Hi(i) = x(i,i)*sum(H_ind);   % Sum row-wose (for each state)
    end
    % Entropy: 
    H = -sum(Hi); 
    return
end

    

% Conditional entropy
% H = sum(p(x,y)log2(p(x)/p(x,y)))
