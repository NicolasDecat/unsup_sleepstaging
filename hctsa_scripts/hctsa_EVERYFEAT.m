
%%%%%%% Together with kmneans_mapping.m and type1auc.m: load all hctsa
%%%%%%% operations, perform kmeans clustering and sequential matching, and
%%%%%%% calculate the one-vs-all type 1 AUC


if isfile('/Users/nico/Documents/HCTSA/Analysis/iterdata.mat') == 1
    delete '/Users/nico/Documents/HCTSA/Analysis/iterdata.mat'
end

tic

Subs = {'807'}; % '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'}; % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};
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

        %% Use feat_id to select data from full op
        datamat = load(hctsafile,'TS_DataMat');
        datamat = datamat.TS_DataMat;
        feat_id = 1:size(datamat,2);  %To include all features

        [timeseries,features]=size(datamat);
        hctsa_ops = datamat(:,feat_id);

        %% Run cross-validation code
        set(0,'DefaultFigureVisible','on') % Remove this to disable the figure displaying (sometimes it could be lots of figures!)
        statistics = [];

        save_stats_columns = {'Type', 'Iteration', 'TrainingAccuracy', 'TestingAccuracy', 'NumberOfFeatures','NumberOfChannels'};
        save_stats = array2table(zeros(0,length(save_stats_columns)));
        save_stats.Properties.VariableNames = save_stats_columns; 
        
        conf = 'BALANCED_LABELED';
        epochSelectFunc = @epochSelect;
        k=5;

        SELECT_TOP_200_FEATURES=size(hctsa_ops,2);
       
        
        [statsOut testMat scoredTest predictTest testTS] = kmeans_mapping_EVERYFEAT(k, hctsa_ops, CM_SAVE_DIR, v, epochSelectFunc, SELECT_TOP_200_FEATURES,sub,v);

        
        scoredTest;  % AASM labels
        predictTest;  % cluster decisions
 
        % Dataset 005: from 381 to 1441
        % Dataset 439: from 151 to 116
        % Dataset 458: from 278 to 1373
        % Dataset 604: from 503 to 1456

        TS = 330:1231;
        scoredTest = scoredTest(TS)';
        predictTest = predictTest(TS)';
        
        Table = table(TS',scoredTest,predictTest,'VariableNames',{'TimeSeries ID','AASM labels','cluster decisions'});



   end

end


