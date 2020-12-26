
%% Computing kmeans clustering and type1 auc for each feature

% Making sure no file remains in folder
if isfile('/Users/nico/Documents/HCTSA/Analysis/AUC/AUC_per_feature.mat') == 1
    delete '/Users/nico/Documents/HCTSA/Analysis/AUC/AUC_per_feature.mat'
end

Subs = {'001'}; % '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
Channels = {'1ch' '2ch' '3ch'};  % used for saveas
NumIter = compose('%diter',(1:100)); % used for saveas


for D = 1:length(Subs)   % For each dataset
    
    sub = Subs{D};

    % File selection
    WHICH_DATA = str2num(sprintf('%s',sub));   

    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))   % get new pwd (current folder)
    HCTSA_DIR = sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub);  % get directory

    ANSWER_FILE=(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_001/ccshs_1800%s_annot.mat',sub));   % annotation file


    
    % Configuration
    addpath '/Users/nico/Documents/GitHub/unsup_sleepstaging';
    configuration_settings;
    clear i n nn op_name name

    
    %% All operation names
    hctsafile = HCTSA_FILE;
    all_op = load(hctsafile,'Operations');
    
    datam = load(hctsafile,'TS_DataMat');   % Load TS_DataMat (epochsxfeatures)

    load('HCTSA_N.mat')
    % for FF = 1:size(all_op.Operations,1)   % For each feature
    
    for v = 1  % For each channel condition
    
        for FF = 1:10
            
            % Load data matrix for one feature
            datamat = datam;
            datamat = datamat.TS_DataMat(:,FF);   % Take data points for 1 feature

            feat_id = 1:size(datamat,2);  

            [timeseries,features]=size(datamat);
            hctsa_ops = datamat(:,feat_id);
            
        
            %% Run cross-validation code
            
            set(0,'DefaultFigureVisible','off') % Remove this to disable the figure displaying 
            exps = EXPS_TO_RUN; 
            statistics = [];

            save_stats_columns = {'Type', 'Iteration', 'TrainingAccuracy', 'TestingAccuracy', 'NumberOfFeatures','NumberOfChannels'};
            save_stats = array2table(zeros(0,length(save_stats_columns)));
            save_stats.Properties.VariableNames = save_stats_columns; 


            if v==1
                NUM_CHANNELS_TO_RUN = [1];
            elseif v == 2
                NUM_CHANNELS_TO_RUN = [2];
            elseif v== 3
                NUM_CHANNELS_TO_RUN = [3];
            end

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

            % Generate all combinations of three channels
            m = [];
            max=200;

            max_test_performance = 0;
            max_combo = [];

            for emg = 10:10:max
                for eog = max-emg:-10:10
                    for eeg = max-emg-eog:-10:10
                        e = [eeg, eog, emg];
                        if sum(e) == 200
                            m = [m; e];
                        end
                    end
                end
            end


            curr_m = [60 110 30];

            SELECT_TOP_200_FEATURES=size(hctsa_ops,2);


            [statsOut testMat scoredTest predictTest Nf testTS AUC_per_feature] = kmeans_mapping_eachfeature(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, SELECT_TOP_200_FEATURES,sub,v,FF,FF);

            [~, statsOut.complexity]=size(hctsa_ops);
            %statsOut.complexity = k;
            statsOut.id = k;
            statistics = [statistics, statsOut];

            iteration=1:size(statsOut.scoredTrain, 1);
            iteration_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.predictTrain)'))/size(statsOut.scoredTrain, 2))';
            iteration_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.predictTest)'))/size(statsOut.scoredTest, 2))';
            iteration_svm_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.svmPredictTrain)'))/size(statsOut.scoredTrain, 2))';
            iteration_svm_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.svmPredictTest)'))/size(statsOut.scoredTest, 2))';
            num_of_features=zeros(size(statsOut.scoredTrain, 1), 1);
            num_of_features(:) = unique(statsOut.totalFeatures);
            num_of_channels=zeros(size(statsOut.scoredTrain, 1), 1);
            num_of_channels(:) = c;


            types=strings(size(statsOut.scoredTrain, 1), 1);
            types(:)=strcat('Unsupervised_', conf);
            row = [types, iteration', iteration_training_accuracy, iteration_testing_accuracy, num_of_features, num_of_channels];
            save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

            types=strings(size(statsOut.scoredTrain, 1), 1);
            types(:)=strcat('Supervised_', conf);
            row = [types, iteration', iteration_svm_training_accuracy, iteration_svm_testing_accuracy, num_of_features, num_of_channels];
            save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

            % end



            %% Draw the confusion matrix for the repeat that has maximum trainCorrect
            if PLOT_CONFUSION_MATRIX
                   [perResponse] = plot_confusion_matrix(statistics.id, statistics.scoredTrain, statistics.predictTrain, ...
                        statistics.scoredTest, statistics.predictTest, CM_SAVE_DIR);
            end

             chan = Channels{NUM_CHANNELS_TO_RUN}; 
    %         fpath = '/Users/nico/Documents/HCTSA/Analysis/CF_v2(10iter)';
    %         saveas(gca,fullfile(fpath,sprintf('CF_%s_%s', sub, chan)),'jpeg')

            % Average percentages of confusion matrices (perResponse)
            % run('Plot_CF_mean')

            set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying
            s=save_stats; [UA, ~, idx] = unique(s(:,[1 6]));NEW_A = [UA,array2table(accumarray(idx,double(table2array(s(:,4))),[],@mean))]; NEW_A;

            % Average percentages of confusion matrices (perResponse)
            Percent_cf(D,v) = {perResponse};   % Row = 1 dataset, col = all iterations 


            %% Plot output accuracy
            accuracy_train = [];
            accuracy_test = [];
            complexity = [];

            for l=1:length(exps)
                k = exps(l);
                accuracy_train = [accuracy_train, statistics(l).output.trainCorrect];
                accuracy_test = [accuracy_test, statistics(l).output.testCorrect];
                complexity = [complexity, statistics(l).complexity];
            end

            FF

          
            
            
        end
     
       
    end

   

    %% Plot average confusion matricees (CF)
    % 
    % % To use for v2 figures
    % Y = cat(3,Percent_cf{:});   % 2nd index = Num channels
    % MEAN_percent_cf = mean(Y,3);  
    % run('Plot_CF_mean.m')
    % 
    % % Average CF of iterations for each channel condition
    % if v<=10   % EEG-only condition
    %     Nit = NumIter{v};    
    %     NUM_CHANNELS_TO_RUN = 1; 
    %     MEAN_percent_cf = mean(Y(:,:,1:v),3);  %  first iterations (CF for EEG only)
    %    run('Plot_CF_mean.m')
    % 
    % elseif v > 10 && v <=20   % EEG+EOG
    %     Nit = NumIter{v/2};  % Half iteration for each of the 2 channels
    %     NUM_CHANNELS_TO_RUN = 1;  % Used for save (chan)
    %     MEAN_percent_cf = mean(Y(:,:,1:10),3);  % Take iterations 1-10 (CF for EEG)
    %     run('Plot_CF_mean.m')
    %     NUM_CHANNELS_TO_RUN = 2; 
    %     MEAN_percent_cf = mean(Y(:,:,11:v),3);  % Take iterations 10-v (CF for EEG+EOG)
    %     run('Plot_CF_mean.m')
    %     
    % elseif v > 20 && v <=30
    % Nit = round(NumIter{v/3});  % third of iteration for each of the 3 channels
    %     NUM_CHANNELS_TO_RUN = 1; 
    %     MEAN_percent_cf = mean(Y(:,:,1:10),3);  % Take iterations 1-10 (CF for EEG)
    %     run('Plot_CF_mean.m')
    %     NUM_CHANNELS_TO_RUN = 2; 
    %     MEAN_percent_cf = mean(Y(:,:,11:20),3);  % Take iterations 11-20 (CF for EEG+EOG)
    %     run('Plot_CF_mean.m')
    %     NUM_CHANNELS_TO_RUN = 3; 
    %     MEAN_percent_cf = mean(Y(:,:,21:v),3);  % Take iterations 21-v (CF for EEG+EOG+EMG)
    %     run('Plot_CF_mean.m')

    
    %% Computing AUC for each feature
    
    
    
    
    
    %% Spectral power analysis: Get the ID of epochs that were labeled wrong by the cluster 
   
%      TableEpochs = [];
%      Stages = [0 1 2 3 5];
%      
%      for stg = 1:5
%          
%          idx = find(statsOut.scoredTest == Stages(stg));     % index of orig labels
%          Cidx = find(statsOut.predictTest == Stages(stg));   % index of cluster decisions
%          AGR = idx(find(ismember(idx,Cidx)));                % index of epochs (from TestTS) who were labelled right by cluster
%          DISAGR = idx(find(~ismember(idx,Cidx)));            % index of epochs (from TestTS) who were labelled right by cluster
%          ID_AGR = testTS(AGR);                               % Actual epochs ID labeled right
%          ID_DISAGR = testTS(DISAGR);                         % Actual epochs ID labeled wrong
%          
%          TableEpochs{stg,1} = idx;   
%          TableEpochs{stg,2} = Cidx;  
%          TableEpochs{stg,3} = AGR;        
%          TableEpochs{stg,4} = DISAGR;       
%          TableEpochs{stg,5} = ID_AGR;                           
%          TableEpochs{stg,6} = ID_DISAGR;  
%                       
%      end
% 
%      % Get the signal of all (dis)agreed epochs. For 2 / 3 channels, index
%      % for TimeSeries.Data is incremented to reach rows for EOG and EMG
%      
%      if v == 1
%         Data = TimeSeries.Data(1:size(eeg_ops,1),:);
%     elseif v == 2
%         Data = TimeSeries.Data((size(eeg_ops,1)+1):(size(eeg_ops,1)*2),:);
%     elseif v == 3
%         Data = TimeSeries.Data((((size(eeg_ops,1)*2)+1):(size(eeg_ops,1)*3)),:);
%      end
%      
%      % 0vs1
%      % 0vs2
%      % 0vs3
%      % 0vs5
%      % 1vs2
%      % 1vs3
%      % 1vs5
%      % 2vs3
%      % 2vs5
%      % 3vs5
% 
% %     % Signal of epochs labelled right / wrong
% %     Signalagr = Data(EpochIDagr,:);
% %     Signaldisagr = TimeSeries.Data(EpochIDdisagr,:);



end


% SummaryTable = table(Iteration,NumChannels,Dataset,Sleep_stage,Testing_accuracy,AUC);
% save('Summary_table_v3.mat','SummaryTable')
% save('CF_data','Percent_cf')

