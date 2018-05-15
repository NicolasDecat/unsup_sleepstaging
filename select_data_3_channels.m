% Select features from reduced_ops.txt, make another .mat file to work on
% cross-validation.
configuration_settings;

homedir = pwd;
%% Start HCTSA tools
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
hctsa_ops = datamat(:,feat_id);

%% Run cross-validation code
% Change the number of operations
set(0,'DefaultFigureVisible','off') % Remove this to disable the figure displaying (sometimes it could be lots of figures!)
exps = EXPS_TO_RUN; % This allow us to selectively choose which experiment to run
channel_statistics = [];
algo_performance  = cell2table(cell(0,6), 'VariableNames', {'dataset', 'run', 'channel', 'exp', 'training_percentage', 'testing_percentage'});

for repeat = 1:10
for channels_used = 1:3 % channels
statistics = [];

for l  = 1:length(exps) 
    k = exps(l); % k is the condition to select operation
    if k==0
        hctsa_ops = datamat(:,feat_id(1:1));
    elseif k==1
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
    statsOut = cross_validation(k, hctsa_ops, CM_SAVE_DIR, channels_used);
    [~, statsOut.complexity]=size(hctsa_ops);
    %statsOut.complexity = k;
    statsOut.id = k;
    statistics = [statistics, statsOut];
    
    row = {WHICH_DATA, repeat, channels_used, k, statsOut.output.trainCorrect, statsOut.output.testCorrect};
    algo_performance = [algo_performance; row];
end

channel_statistics = [channel_statistics, statistics];
    
end
end

save('ALGO_PERFORMANCE.mat', 'algo_performance');

%% Draw the confusion matrix for the repeat that has maximum trainCorrect
% if PLOT_CONFUSION_MATRIX
%     for idx = 1:length(exps)
%         plot_confusion_matrix(statistics(idx).id, statistics(idx).scoredTrain, statistics(idx).predictTrain, ...
%             statistics(idx).scoredTest, statistics(idx).predictTest, CM_SAVE_DIR)
%     end
% end
% 
% set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying
% 
% figure;
% 
% for c = 1:3 % channels
% %% Plot output accuracy
% accuracy_train = [];
% accuracy_test = [];
% complexity = [];
% 
% for l=1:length(exps)
%     k = exps(l);
%     stats=channel_statistics((length(exps)*(c-1))+l);
%     accuracy_train = [accuracy_train, stats.output.trainCorrect];
%     accuracy_test = [accuracy_test, stats.output.testCorrect];
%     complexity = [complexity, stats.complexity];
% end
% 
% semilogx(complexity,accuracy_test);
% hold on;
% 
% end
% 
% title(sprintf('Accuracy Report between channels for dataset %d', WHICH_DATA));
% ylabel('Accuracy [0-1]');
% xlabel('Number of features');
% xlim([1 200]);
% 
% legend('1 channel','2 channels', '3 channels');
% hold off;

