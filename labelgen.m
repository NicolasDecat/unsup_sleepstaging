function [labels,keywords] = labelgen(no_timeseries,no_level,labelname1,labelname2,labelname3);
% ===============================================================
% Function to automatically label time series data according to the input
% label name, 3 levels of subgrouping
% Created: 20/12/16
% Modified: 9/1/17
% ###############################################################
% Input arguments: no_timeseries - number of time series in data matrix
%                  no_level - number of groups/subgroups
%                  labelname - cell array of label for each group/subgroup
% Output: Labels - cell array of label for each time series
%         Keywords - cell array of keywords for each time series
% ################################################################
% Check input arguments 
if nargin<1
    error('Insufficient input. Must enter number of time series.')
end

if nargin<2 || isempty(no_level)
    error('Must enter number of subgroup levels.')
end

if nargin<3 || isempty(labelname1)
    error('Insufficient inputs, must enter labels of first group.')
end

if nargin<4 && no_level<2
   labelname2=[];
   labelname3=[];
elseif (nargin<4 || isempty(labelname2)) && no_level>=2
    error('Insufficient input, must enter labels of second group.')
end

if nargin<5 && no_level<3
   labelname3=[];
elseif (nargin<5 || isempty(labelname3)) && no_level>=3
    error('Insufficient input, must enter labels of third group.')
end   


% Check if labelname are cell array and not empty (later) 
if ~iscellstr(labelname1)
    error('Subgroup 1 labels are invalid.')
elseif ~iscellstr(labelname2) && no_level>1
    error('Subgroup 2 labels are invalid.')
elseif ~iscellstr(labelname3) && no_level>2
    error('Subgroup 3 labels are invalid.')
end
% Number of groups
n_group1 = length(labelname1);
n_group2 = length(labelname2);
n_group3 = length(labelname3);

% Number of time series per group
n_ts1 = no_timeseries/n_group1;
n_ts2 = n_ts1/n_group2;
n_ts3 = n_ts2/n_group3;


% Assign labels and keywords to each time series
% without checking conditions, assume all data provided

switch no_level
    case 1
        for i=1:n_group1
            labels(n_ts1*(i-1)+1:n_ts1*i)= labelname1(i);
            keywords(n_ts1*(i-1)+1:n_ts1*i) = labelname1(i);
        end
    case 2
        for i=1:n_group1
            labels(n_ts1*(i-1)+1:n_ts1*i)= labelname1(i);
            keywords(n_ts1*(i-1)+1:n_ts1*i) = labelname1(i);
            star = n_ts1*(i-1);
            for j=1:n_group2
                labels(n_ts2*(j-1)+1+star:n_ts2*j+star)= strcat(labels(n_ts2*(j-1)+1+star:n_ts2*j+star),'_',repmat(labelname2(j),size(labels(n_ts2*(j-1)+1+star:n_ts2*j+star))));
                keywords(n_ts2*(j-1)+1+star:n_ts2*j+star) = strcat(keywords(n_ts2*(j-1)+1+star:n_ts2*j+star),',',repmat(labelname2(j),size(labels(n_ts2*(j-1)+1+star:n_ts2*j+star))));
            end 
        end
    case 3
        for i=1:n_group1
            labels(n_ts1*(i-1)+1:n_ts1*i)= labelname1(i);
            keywords(n_ts1*(i-1)+1:n_ts1*i) = labelname1(i);
            star = n_ts1*(i-1);
            for j=1:n_group2
                labels(n_ts2*(j-1)+1+star:n_ts2*j+star)= strcat(labels(n_ts2*(j-1)+1+star:n_ts2*j+star),'_',repmat(labelname2(j),size(labels(n_ts2*(j-1)+1+star:n_ts2*j+star))));
                keywords(n_ts2*(j-1)+1+star:n_ts2*j+star) = strcat(keywords(n_ts2*(j-1)+1+star:n_ts2*j+star),',',repmat(labelname2(j),size(labels(n_ts2*(j-1)+1+star:n_ts2*j+star))));
                moon =n_ts2*(j-1)+star;
                for k=1:n_group3
                    labels(n_ts3*(k-1)+1+moon:n_ts3*k+moon)= strcat(labels(n_ts3*(k-1)+1+moon:n_ts3*k+moon),'_',repmat(labelname3(k),size(labels(n_ts3*(k-1)+1+moon:n_ts3*k+moon))));
                    keywords(n_ts3*(k-1)+1+moon:n_ts3*k+moon) = strcat(keywords(n_ts3*(k-1)+1+moon:n_ts3*k+moon),',',repmat(labelname3(k),size(labels(n_ts3*(k-1)+1+moon:n_ts3*k+moon))));
                end
            end 
        end
end
return