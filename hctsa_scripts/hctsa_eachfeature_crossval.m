
%%%%%%%% Here is the correct verion: Together with
%%%%%%%% kmeans_mapping_eachfeature_crossval.m and
%%%%%%%% type1auc_eachfeature_crossval.m: we compute accuracy for each
%%%%%%%% feature (using leave-1-out strategy and threshold) and each binary classifier instead of using kmeans and AUC. 


%% Computing kmeans clustering and type1 auc for each feature

tic

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
Channels = {'1ch' '2ch' '3ch'};  % used for saveas
NumIter = compose('%diter',(1:100)); % used for saveas


for D = 1:length(Subs)   % For each dataset
    
    % Making sure no file remains in folder
    if isfile('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Per_correct_mean.mat') == 1
        delete '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Per_correct_mean.mat'
    end
    
    sub = Subs{D};

    % Configuration
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))   % get new pwd (current folder)
    addpath '/Users/nico/Documents/GitHub/unsup_sleepstaging';
    configuration_settings;
    clear i n nn op_name name

    
    %% Load all operations
    hctsafile = HCTSA_FILE;
    all_op = load(hctsafile,'Operations');
    datam = load(hctsafile,'TS_DataMat');   % Load TS_DataMat (epochsxfeatures)
    load('HCTSA_N.mat')
    
    
    % Take 5%
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Top_5perc')
    idx_feat = cell2mat((Top_5perc(:,1)));

    % Equivalent index 
    for i=1:numel(idx_feat)
        BestFeat(i) = find(Operations.ID == idx_feat(i));
    end
               
    
    for v = 1:1  % For each channel condition

        %for FF = 1:size(Operations,1)    % For each feature
        for FF = 1:280
       
          
            % Load data matrix for one feature
            datamat = datam;
            % datamat = datamat.TS_DataMat(:,FF);   % Take data points for 1 feature
            datamat = datamat.TS_DataMat(:,BestFeat(FF));   % Take data points for 1 feature

            feat_id = 1:size(datamat,2);  

            [timeseries,features]=size(datamat);
            hctsa_ops = datamat(:,feat_id);
          
        
            %% Run cross-validation code
            
            set(0,'DefaultFigureVisible','on') % Remove this to disable the figure displaying 
            exps = EXPS_TO_RUN; 
            statistics = [];

            save_stats_columns = {'Type', 'Iteration', 'TrainingAccuracy', 'TestingAccuracy', 'NumberOfFeatures','NumberOfChannels'};
            save_stats = array2table(zeros(0,length(save_stats_columns)));
            save_stats.Properties.VariableNames = save_stats_columns; 

            % Channel(s) to use
            NUM_CHANNELS_TO_RUN = v;
            chan = Channels{NUM_CHANNELS_TO_RUN};  % used for saveas
            c = NUM_CHANNELS_TO_RUN;
            
            l = 1:length(exps);
            k = exps(l); 
            hctsa_ops = datamat(:,feat_id);

            conf = 'BALANCED_LABELED';
            epochSelectFunc = @epochSelect;


            single_channel_size=size(hctsa_ops,1)/NUM_CHANNELS;
            eeg_ops=hctsa_ops(1:single_channel_size,:);
            eog_ops=hctsa_ops(single_channel_size+1:single_channel_size*2,:);
            emg_ops=hctsa_ops(single_channel_size*2+1:single_channel_size*3,:);


            SELECT_TOP_200_FEATURES=size(hctsa_ops,2);


            [testMat Nf testTS Per_correct_mean] = kmeans_mapping_eachfeature_crossval(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, SELECT_TOP_200_FEATURES,sub,v,FF,FF);

            % [~, statsOut.complexity]=size(hctsa_ops);
            
            %%%%% statsOut.complexity = k;
            
            statsOut.id = k;
            statistics = [statistics, statsOut];

%             iteration=1:size(statsOut.scoredTrain, 1);
%             iteration_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.predictTrain)'))/size(statsOut.scoredTrain, 2))';
%             iteration_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.predictTest)'))/size(statsOut.scoredTest, 2))';
% %%%%%             iteration_svm_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.svmPredictTrain)'))/size(statsOut.scoredTrain, 2))';
% %%%%%             iteration_svm_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.svmPredictTest)'))/size(statsOut.scoredTest, 2))';
%             num_of_features=zeros(size(statsOut.scoredTrain, 1), 1);
%             num_of_features(:) = unique(statsOut.totalFeatures);
%             num_of_channels=zeros(size(statsOut.scoredTrain, 1), 1);
%             num_of_channels(:) = c;


%             types=strings(size(statsOut.scoredTrain, 1), 1);
%             types(:)=strcat('Unsupervised_', conf);
%             row = [types, iteration', iteration_training_accuracy, iteration_testing_accuracy, num_of_features, num_of_channels];
%             save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];
% 
%             types=strings(size(statsOut.scoredTrain, 1), 1);
%             types(:)=strcat('Supervised_', conf);
%             row = [types, iteration', iteration_svm_training_accuracy, iteration_svm_testing_accuracy, num_of_features, num_of_channels];
%             save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

            % end



%             %% Draw the confusion matrix for the repeat that has maximum trainCorrect
%             if PLOT_CONFUSION_MATRIX
%                    [perResponse] = plot_confusion_matrix(statistics.id, statistics.scoredTrain, statistics.predictTrain, ...
%                         statistics.scoredTest, statistics.predictTest, CM_SAVE_DIR);
%             end

             %chan = Channels{NUM_CHANNELS_TO_RUN}; 
    %         fpath = '/Users/nico/Documents/HCTSA/Analysis/CF_v2(10iter)';
    %         saveas(gca,fullfile(fpath,sprintf('CF_%s_%s', sub, chan)),'jpeg')

            % Average percentages of confusion matrices (perResponse)
            % run('Plot_CF_mean')

             %set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying
             % s=save_stats; [UA, ~, idx] = unique(s(:,[1 6]));NEW_A = [UA,array2table(accumarray(idx,double(table2array(s(:,4))),[],@mean))]; NEW_A;

            % Average percentages of confusion matrices (perResponse)
%             Percent_cf(D,v) = {perResponse};   % Row = 1 dataset, col = all iterations 


            %% Plot output accuracy
%             accuracy_train = [];
%             accuracy_test = [];
%             complexity = [];
% 
%             for l=1:length(exps)
%                 k = exps(l);
%                 accuracy_train = [accuracy_train, statistics(l).output.trainCorrect];
%                 accuracy_test = [accuracy_test, statistics(l).output.testCorrect];
%                 complexity = [complexity, statistics(l).complexity];
%             end


        end
   
       
    end
    
    gpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_Top_5perc/average_univariate(280)';
    save(fullfile(gpath,sprintf('Per_correct_mean(Dataset %s).mat',sub)),'Per_correct_mean')  % Save all columns in AUC folder


    %% Plot Data matrix

%     figure;
%     imagesc(Per_correct_mean)
%     title(sprintf('Classification performance per feature (Dataset %s)',sub));
%     % title('Classification performance per feature (Dataset 001)');
% 
%     ax = gca;
%     ax.XTick = 1:500:size(Operations,1);
%     ax.YTick = 1:10;
%     % ax.XTickLabels = strseq('f',1:100)';
%     ax.XTickLabels = arrayfun(@(a)num2str(a),0:500:size(Operations,1),'uni',0);
%     ax.YTickLabels = {'W vs N1', 'W vs N2', 'W vs N3', 'W vs REM', 'N1 vs N2','N1 vs N3','N1 vs REM','N2 vs N3','N2 vs REM','N3 vs REM'};
%     % ax.YTickLabels = {'Wake vs all', 'N1 vs all', 'N2 vs all', 'N3 vs all', 'REM vs all'};
%     ylabel('Binary classifiers');
%     xlabel('features');
%     ax.XAxisLocation = 'bottom';
% 
%     colormap 'default'
%     colorbar
% 
% %    % Save figure
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Figure_accuracy_per_feat';
% %     saveas(gca,fullfile(fpath,sprintf('AccperFeat(100)_%s_crossval',sub)),'fig')
% %     saveas(gca,fullfile(fpath,sprintf('AccperFeat(100)_%s_crossval',sub)),'jpg')
%  
% %   % Save matrix
% %     gpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat';
% %     save(fullfile(gpath,sprintf('Per_correct_mean_OVA(Dataset %s).mat',sub)),'Per_correct_mean') 
%        gpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy/Matrix_accuracy_per_feat/OVA';
% %        save(fullfile(gpath,sprintf('Per_correct_mean_OVA(Dataset %s).mat',sub)),'Per_correct_mean')


end

toc

