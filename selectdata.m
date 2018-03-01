% Select features from reduced_ops.txt, make another .mat file to work on
% cross-validation.
clear all; clc;

configuration_settings;

homedir = pwd;
% Start HCTSA tools
cd(HCTSA_DIR)
startup
cd(homedir)

%% Read text file
fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);

%% Wanted operation names
feat_name = features{1,2};

%% All operation names
hctsafile = HCTSA_FILE;
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

% annotation = load(ANSWER_FILE);
% label = annotation.sleepstage;
% 
%% Perform feature selection (experimental)

% Include dependencies
% addpath(strcat(FSLIB_TOOLBOX_DIR, filesep, 'lib')); % dependencies
% addpath(strcat(FSLIB_TOOLBOX_DIR, filesep, 'methods')); % FS methods
% addpath(genpath(strcat(FSLIB_TOOLBOX_DIR, filesep, 'lib/drtoolbox')));
% 
% features_channel1 = datamat(334:1374, :);
% features_channel2 = datamat((1374+334):1374*2, :);
% features_channel3 = datamat(((1374*2)+334):1374*3, :);

% X_train = [features_channel1; features_channel2; features_channel3];
% Y_train = [label(334:1374); label(334:1374); label(334:1374)];

% X_train = features_channel1;
% Y_train = label(334:1374);
% numF = size(X_train, 2);

%[ranking, w] = reliefF(x_data, y_data, 20);
%[ranking, w, subset] = ILFS_auto(x_data, y_data, 4, 0 )
%ranking = mRMR(x_data, y_data, size(x_data, 2));
%[ ranking , w] = mutInfFS( X_train, Y_train, size(X_train , 2));
%[ ranking , w] = fsvFS( X_train, Y_train, size(X_train , 2) );

%Laplacian
% W = dist(X_train');
% W = -W./max(max(W)); % it's a similarity
% [lscores] = LaplacianScore(X_train, W);
% [junk, ranking] = sort(-lscocd cd gires);

% MCFS: Unsupervised Feature Selection for Multi-Cluster Data
% options = [];
% options.k = 5; %For unsupervised feature selection, you should tune
% %this parameter k, the default k is 5.
% options.nUseEigenfunction = 4;  %You should tune this parameter.
% [FeaIndex,~] = MCFS_p(X_train,numF,options);
% ranking = FeaIndex{1};
        
% ranking = spider_wrapper(X_train,Y_train,numF,lower('rfe'));
% ranking = spider_wrapper(X_train,Y_train,numF,lower('10'));
% ranking = spider_wrapper(X_train,Y_train,numF,lower('fisher'));

% Infinite Feature Selection 2015 updated 2016
% alpha = 0.5;    % default, it should be cross-validated.
% sup = 1;        % Supervised or Not
% [ranking, w] = infFS( X_train , Y_train, alpha , sup , 0 );    

% This is matlab feature selection
%[ranked, weight] = relieff(features_channel1, label(trainTS), 10, 'method', 'classification', 'categoricalx', 'on');

% Features Selection via Eigenvector Centrality 2016
% alpha = 0.5; % default, it should be cross-validated.
% ranking = ECFS( X_train, Y_train, alpha )  ;

% Regularized Discriminative Feature Selection for Unsupervised Learning
% nClass = 2;
% ranking = UDFS(X_train , nClass ); 

% BASELINE - Sort features according to pairwise correlations
% ranking = cfs(X_train);     
%         
% [B, I] = sort(ranking);
% feat_id = I;


%% Run cross-validation code
% Change the number of operations
%set(0,'DefaultFigureVisible','off') % Remove this to disable the figure displaying (sometimes it could be lots of figures!)
for k = 1:10 % k is the condition to select operation
%for k = 5:5 % k is the condition to select operation
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
        %hctsa_ops = datamat(:,feat_id(1:200));
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
%set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying

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

saveas(gcf, strcat(CM_SAVE_DIR, filesep, 'ACCURACY_REPORT.png'));

%% Checking 

% for n=1:100
%     a = epochSampling(14);
%     boundary(n,:) = [min(a), max(a)];
% end
% minmin = min(boundary(:,1))
% maxmax = max(boundary(:,2))

