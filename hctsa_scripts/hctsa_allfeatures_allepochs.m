
%%%%%%% Together with kmneans_mapping.m and type1auc.m: load all hctsa
%%%%%%% operations, perform kmeans clustering and sequential matching, and
%%%%%%% calculate the one-vs-all type 1 AUC


if isfile('/Users/nico/Documents/HCTSA/Analysis/iterdata.mat') == 1
    delete '/Users/nico/Documents/HCTSA/Analysis/iterdata.mat'
end

tic

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
Channels = {'1ch' '2ch' '3ch'};  % used for save
NumIter = compose('%diter',(1:100)); % used for save


for D = 1:length(Subs)
    
    sub = Subs{D};

    % File selection
    WHICH_DATA = str2num(sprintf('%s',sub)); 

    % go to current folder
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    HCTSA_DIR = sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub);  
    ANSWER_FILE=(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s/ccshs_1800%s_annot.mat',sub,sub));

    
    for v = 3   % For each channel condition

        % Configuration
        addpath '/Users/nico/Documents/GitHub/unsup_sleepstaging';
        configuration_settings;

        %% All operation names
        hctsafile = HCTSA_FILE;
        all_op = load(hctsafile,'Operations');

        clear i n nn op_name name

        %% Use feat_id to select data from full op
        datamat = load(hctsafile,'TS_DataMat');
        datamat = datamat.TS_DataMat;
        feat_id = 1:size(datamat,2);  %To include all features

        [timeseries,features]=size(datamat);
        hctsa_ops = datamat(:,feat_id);

        %% Run cross-validation code
        % Change the number of operations
        set(0,'DefaultFigureVisible','off') % Remove this to disable the figure displaying (sometimes it could be lots of figures!)
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

        chan = Channels{NUM_CHANNELS_TO_RUN};  % used for save

        c = NUM_CHANNELS_TO_RUN;
        l = 1:length(exps);
        k = exps(l); % k is the condition to select operation
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
       
        
        [statsOut testMat scoredTest predictTest Nf Iteration NumChannels Dataset Sleep_stage Testing_accuracy AUC testTS] = kmeans_mapping_allepochs(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, SELECT_TOP_200_FEATURES,sub,v);
        [~, statsOut.complexity]=size(hctsa_ops);
        %statsOut.complexity = k;
        statsOut.id = k;
        statistics = [statistics, statsOut];

        iteration=1:size(statsOut.scoredTrain, 1);
        iteration_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.predictTrain)'))/size(statsOut.scoredTrain, 2))';
        iteration_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.predictTest)'))/size(statsOut.scoredTest, 2))';
%         iteration_svm_training_accuracy = ((sum((statsOut.scoredTrain == statsOut.svmPredictTrain)'))/size(statsOut.scoredTrain, 2))';
%         iteration_svm_testing_accuracy = ((sum((statsOut.scoredTest == statsOut.svmPredictTest)'))/size(statsOut.scoredTest, 2))';
        num_of_features=zeros(size(statsOut.scoredTrain, 1), 1);
        num_of_features(:) = unique(statsOut.totalFeatures);
        num_of_channels=zeros(size(statsOut.scoredTrain, 1), 1);
        num_of_channels(:) = c;


        types=strings(size(statsOut.scoredTrain, 1), 1);
        types(:)=strcat('Unsupervised_', conf);
        row = [types, iteration', iteration_training_accuracy, iteration_testing_accuracy, num_of_features, num_of_channels];
        save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

%         types=strings(size(statsOut.scoredTrain, 1), 1);
%         types(:)=strcat('Supervised_', conf);
%         row = [types, iteration', iteration_svm_training_accuracy, iteration_svm_testing_accuracy, num_of_features, num_of_channels];
%         save_stats = [save_stats; array2table(row, 'VariableNames', save_stats_columns)];

        % end



        %% Draw the confusion matrix for the repeat that has maximum trainCorrect
        if PLOT_CONFUSION_MATRIX
               [perResponse] = plot_confusion_matrix(statistics.id, statistics.scoredTrain, statistics.predictTrain, ...
                    statistics.scoredTest, statistics.predictTest, CM_SAVE_DIR);
        end
    
        chan = Channels{NUM_CHANNELS_TO_RUN}; 
        fpath = '/Users/nico/Documents/HCTSA/Analysis/AUC_100/ConfusionMatrix';
        saveas(gca,fullfile(fpath,sprintf('CF_%s_%s', sub, chan)),'jpeg')
        
        set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying
        s=save_stats; [UA, ~, idx] = unique(s(:,[1 6]));NEW_A = [UA,array2table(accumarray(idx,double(table2array(s(:,4))),[],@mean))]; NEW_A;
        
        % Average percentages of confusion matrices (perResponse)
        Percent_cf(D,v) = {perResponse};   % Row = 1 dataset, col = channel (mean over iterations) 
        

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
  

    end
   
    fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms/allepochs';
    save([fpath sprintf('/statsOut_allepochs_3ch(%s)',sub)],'statsOut')

end

%%% Average CF (already averaged over iterations) over datasets for each channel condition

% EEG only
% CH = 1;
% Y = cat(3,Percent_cf{:,1});   
% MEAN_percent_cf = mean(Y,3);
% run('Plot_CF_mean.m') % change saveas name
% 
% % EEG+EOG
% CH = 2;
% Y = cat(3,Percent_cf{:,2});   
% MEAN_percent_cf = mean(Y,3);
% run('Plot_CF_mean.m')
% 
% % EEG+EOG+EMG
% CH = 3;
% Y = cat(3,Percent_cf{:,3});   
% MEAN_percent_cf = mean(Y,3);
% run('Plot_CF_mean.m')
% 
% % All
% Y = cat(3,Percent_cf{:});   
% MEAN_percent_cf = mean(Y,3);
% run('Plot_CF_mean.m')


SummaryTable = table(Iteration,NumChannels,Dataset,Sleep_stage,Testing_accuracy,AUC);
% save('Summary_table_v3.mat','SummaryTable')
% save('CF_data','Percent_cf')

toc

%%%% Plot TS closest to k centroid
% 
% for C = 1:5   % for each cluster
%    
%     [~,I] = sort(D(:,C),'ascend');   % sort, from closest to furthest from centroid
%     closest(C) = I(1);               % get the closest TS
%     
%     %plot
%     TS = TimeSeries.Data{closest(C),:};
%     figure;plot(TS)
%     
% end
% 

