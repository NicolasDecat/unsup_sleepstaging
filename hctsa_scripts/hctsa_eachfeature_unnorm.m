

%% Computing kmeans clustering and type1 auc for each feature

tic

% Subs = {'439'};  % '001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821'};
Subs = {'001' '005' '458' '596' '748' '749' '752' '604' '807' '821'};

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
    
    % Unnormalised features, 5603 features
    load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/ID5603')
    
    load('HCTSA.mat')
    Operations = Operations(ID5603,:);
    TS_DataMat = TS_DataMat(:,ID5603);
    
             
    for v = 1:1  % For each channel condition

        for FF = 1:size(Operations,1)    % For each feature
       
            % Load data matrix for one feature
            datamat = datam;
            datamat = datamat.TS_DataMat(:,FF);   % Take data points for 1 feature

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


            SELECT_TOP_200_FEATURES=size(hctsa_ops,2);


            [testMat Nf testTS Per_correct_mean] = kmeans_mapping_eachfeature_unnorm(k, hctsa_ops, CM_SAVE_DIR, c, epochSelectFunc, SELECT_TOP_200_FEATURES,sub,v,FF,FF);

            statsOut.id = k;
            statistics = [statistics, statsOut];


        end
   
       
    end
    
    gpath = '/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_each_unnorm';
    save(fullfile(gpath,sprintf('Per_correct_mean(Dataset %s).mat',sub)),'Per_correct_mean')  % Save all columns in AUC folder


end


%% Mean rank of ind features

load('HCTSA.mat', 'Operations')

% All datasets
Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821'};
    
for D = 1:11

    sub = Subs{D};

    % Rank features and get their ID, names, accuracy
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/unsup_each_unnorm/Per_correct_mean(Dataset %s)',sub))

    Matrix{D} = Per_correct_mean;
    
end

Y = cat(3,Matrix{:});
Matrix_allD = mean(Y,3);

means = mean(Matrix_allD);
[~,I] = sort(means','descend');
means = means(I);

Matrix_allD = Matrix_allD(:,I);
Top_Feat = I(1:5603); 

% Get the name and keyword associated with these features
CodeString = {Operations.CodeString}.';
Keywords = {Operations.Keywords}.';
YLabel = {Operations.ID}.';

Top_name(1:5603) = CodeString(Top_Feat,1);
Top_key(1:5603) = Keywords(Top_Feat,1);
Top_ID(1:5603) = YLabel(Top_Feat,1);

Top_5603 = [Top_name' Top_key' Top_ID'];
% 
%     % Top_mean (mean over classifiers)
%     for F = 1:length(Top_5603)
%         for C = 1:10
%             Accu(F,C) = mean(Per_correct_mean(C,F));
%         end
%     end
% 
%     Top_mean = mean(Accu');
%         
%     clear top_5603
    

 