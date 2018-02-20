% Select features from reduced_ops.txt, make another .mat file to work on
% cross-validation.
clear all; clc;
%% Read text file
fileID = fopen('reduced_ops.txt');
features = textscan(fileID,'%s %s %s');
fclose(fileID);

%% Wanted operation names
feat_name = features{1,2};

%% All operation names
hctsafile = 'HCTSA_N.mat';
all_op = load(hctsafile,'Operations');


%% Check operation name, get feat_id
nn=0;
for n = 1:length(feat_name)
    op_name = char(feat_name(n));
    for i = 1:length(all_op.Operations)
        name = all_op.Operations(i).Name;
        if strcmp(op_name,name)
            nn=nn+1;
            feat_id(nn) = i;
            feat(nn).id = i; % all_op.Operations(i).ID % Actual operation ID
            feat(nn).name = cellstr(name);
        end
    end
end
clear i n nn op_name name
%% Use feat_id to select data from full op
datamat = load(hctsafile,'TS_DataMat');
datamat = datamat.TS_DataMat;

[timeseries,features]=size(datamat);
% hctsa_ops = datamat(:,feat_id);
% v
% v
% v
%% Run cross-validation code
% Change the number of operations
for k = 1:10 % k is the condition to select operation
    if k==1
        hctsa_ops = datamat(:,feat_id(1:10));
    elseif k==2
        hctsa_ops = datamat(:,feat_id(1:25));
    elseif k==3
        hctsa_ops = datamat(:,feat_id(1:50));
    elseif k==4
        hctsa_ops = datamat(:,feat_id(1:100));
    elseif k==5 % Top 198 features
        hctsa_ops = datamat(:,feat_id);
    elseif k==6 % Random 500 features
        rand_id = randperm(features,500);
        hctsa_ops = datamat(:,rand_id);
    elseif k==7
        rand_id = randperm(features,1000);
        hctsa_ops = datamat(:,rand_id);
    elseif k==8
        rand_id = randperm(features,2000);
        hctsa_ops = datamat(:,rand_id);
    elseif k==9
        rand_id = randperm(features,5000);
        hctsa_ops = datamat(:,rand_id);
    else
        hctsa_ops = datamat;
    end
    % run('crossvalKR.m') 
    [~,complexity(k)]=size(hctsa_ops);
    run crossval.m
end

%% Run 1 condition
hctsa_ops = datamat(:,feat_id);
k=1;
%run('crossval.m')

%% Plot output accuracy
for k=1:length(Output)
    accuracy_train(k) = Output(k).trainCorrect;
    accuracy_test(k) = Output(k).testCorrect;
end

figure;
semilogx(complexity,accuracy_train,complexity,accuracy_test)
legend('Training','Test')
ylabel('Accuracy [0-1]')
xlabel('Number of features')


%% Checking 

for n=1:100
    a = epochSampling(14);
    boundary(n,:) = [min(a), max(a)];
end
minmin = min(boundary(:,1))
maxmax = max(boundary(:,2))

