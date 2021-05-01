
%%%%%%%% Goes with the crossval script group: computing the accuracy using
%%%%%%%% all features at once for each binary classifier (kmeans returning
%%%%%%%% 2 clusters for each classifier), to compare the difference in
%%%%%%%% accuracy when using all features at once and only 1 feature at a
%%%%%%%% time (the latter being returned by crossval script group)


%% Computing kmeans clustering and type1 auc for each feature

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
Channels = {'1ch' '2ch' '3ch'};  % used for saveas
NumIter = compose('%diter',(1:100)); % used for saveas

 
for D = 1:length(Subs)
    
    sub = Subs{D};

    % File selection
    WHICH_DATA = str2num(sprintf('%s',sub)); 

    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) % get new pwd (current folder)
    HCTSA_DIR = sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub);  % get directory

    ANSWER_FILE=(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s/ccshs_1800%s_annot.mat',sub,sub));

    
    for v = 1:1   % For each channel condition

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
        
        Catch22_idx = [11,12,134,135,243,1121,7634,3477,1406,1585,1965,310,2997,3264,3294,4492,3467,3604,4036,4156,4421,3010];
       
        % Convert 7749->6006
        ID = all_op.Operations{:,4};
        
        clear Equi_Catch22
        Equi_Catch22 = find(ismember(ID,Catch22_idx));

        feat_id = Equi_Catch22;  %To include all features

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


        [testMat Nf testTS trainTS trainMat label stats PERC_STAGE] = kmeans_catch22(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, SELECT_TOP_200_FEATURES,sub,v,D,D);

        
%          % SUPERVISED (same data)
%         template = templateSVM(...
%         'KernelFunction', 'linear', ...
%         'PolynomialOrder', [], ...
%         'KernelScale', 'auto', ...
%         'BoxConstraint', 1, ...
%         'Standardize', true);
% 
%     
%          for T = 1:length(trainMat)
%              
%              SVMModel{T} = fitcecoc(trainMat{T}, label(trainTS),'Learners', template, ...
%              'Coding', 'onevsone', 'ClassNames', [0; 1; 2; 3; 5]);
%      
%              svmTrain(:,T) = predict(SVMModel{T}, trainMat{T});
%              svmTest(:,T) = predict(SVMModel{T}, testMat{T});
%              
%              stats.svmPredictTrain(T, :) = svmTrain(:,T)';  
%              stats.svmPredictTest(T, :) = svmTest(:,T)';
%              
%              iteration_svm_training_accuracy(T,:) = ((sum((stats.scoredTrain(T,:) == stats.svmPredictTrain(T,:))'))/size(stats.scoredTrain, 2))';
%              iteration_svm_testing_accuracy(T,:) = ((sum((stats.scoredTest(T,:) == stats.svmPredictTest(T,:))'))/size(stats.scoredTest, 2))';
%          end
%          
%          iteration_svm_training_accuracy_mean = mean(iteration_svm_training_accuracy);
%          iteration_svm_testing_accuracy_mean = mean(iteration_svm_testing_accuracy);
% 
%   
    end
    
    
    PERC_PER_CLASSIFIER_10iter = mean(PERC_STAGE); 
    
    % Save for each dataset
    fpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_catch22';
    save(fullfile(fpath,sprintf('PERC_PER_CLASSIF_10iter_Dataset%s.mat',sub)),'PERC_PER_CLASSIFIER_10iter')  

    
end


