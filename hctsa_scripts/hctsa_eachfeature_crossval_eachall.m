
%%%%%%%% Goes with the crossval script group: computing the accuracy using
%%%%%%%% all features at once for each binary classifier (kmeans returning
%%%%%%%% 2 clusters for each classifier), to compare the difference in
%%%%%%%% accuracy when using all features at once and only 1 feature at a
%%%%%%%% time (the latter being returned by crossval script group)


%% Computing kmeans clustering and type1 auc for each feature


Subs = {'001'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
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
        feat_id = 1:size(datamat,2);  %To include all features

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


        [testMat Nf testTS PERC_STAGE] = kmeans_mapping_eachfeature_crossval_eachall(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, SELECT_TOP_200_FEATURES,sub,v);

       
    end
    
    
    PERC_PER_CLASSIFIER_10iter = mean(PERC_STAGE); 
    
    % Save for each dataset
    fpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy';
    save(fullfile(fpath,sprintf('PERC_PER_CLASSIFIER_30iter_Dataset%s.mat',sub)),'PERC_PER_CLASSIFIER_10iter')  

    
end


