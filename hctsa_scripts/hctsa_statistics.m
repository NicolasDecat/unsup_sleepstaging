2
%% Compute the mean type-1 AUC

load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');

% Mean AUC for EEG condition
idx_EEG = find(Table_AUC_100.NumChannels == 1);
mean_AUC_EEG = mean(Table_AUC_100.AUC(idx_EEG,1));  % 0.6948
std_AUC_EEG = std(Table_AUC_100.AUC(idx_EEG,1));    % 0.1498

% Mean AUC for EEG+EOG condition
idx_EEGEOG = find(Table_AUC_100.NumChannels == 2);
mean_AUC_EEGEOG = mean(Table_AUC_100.AUC(idx_EEGEOG,1));  % 0.7290
std_AUC_EEGEOG = std(Table_AUC_100.AUC(idx_EEGEOG,1));    % 0.1460

% Mean AUC for EEG+EOG+EMG condition
idx_EEGEOGEMG = find(Table_AUC_100.NumChannels == 3);
mean_AUC_EEGEOGEMG = mean(Table_AUC_100.AUC(idx_EEGEOGEMG,1));  % 0.7433
std_AUC_EEGEOGEMG = std(Table_AUC_100.AUC(idx_EEGEOGEMG,1));    % 0.1334

%% Paired t-test between mean EEG and mean EEG+EOG+EMG (clustering)

%%%% Average across 100 iterations for each dataset (EEG condition)
load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');
idx_EEG = find(Table_AUC_100.NumChannels == 1);
Table_AUC_100 = Table_AUC_100(idx_EEG,:);

Subs = [1 5 439 458 596 748 749 752 604 807 821 870];

% Get mean across iterations for each dataset
for D = 1:12
    a = find(Table_AUC_100.Dataset == Subs(D));
    EEG(D) = mean(Table_AUC_100.AUC(a));    % AUC for EEG condition  
end

%%%% Average across 100 iterations for each dataset (EEG+EOG+EMG condition)
load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');
idx_EEGEOGEMG = find(Table_AUC_100.NumChannels == 3);
Table_AUC_100 = Table_AUC_100(idx_EEGEOGEMG,:);

for D = 1:12
    a = find(Table_AUC_100.Dataset == Subs(D));
    EEGEOGEMG(D) = mean(Table_AUC_100.AUC(a));  % AUC for EEG+EOG+EMG condition
end

% Average and std for both conditions
mean_EEG = mean(EEG);
std_EEG = std(EEG);

mean_EEGEOGEMG = mean(EEGEOGEMG);
std_EEGEOGEMG = std(EEGEOGEMG);

% Paired t-test
[h,p,ci,stats] = ttest(EEGEOGEMG,EEG)

% Paired t-test NOV
[h,p,ci,stats] = ttest(EEGEOGEMG,EEG)


% figure; plot(EEG)
% hold on 
% plot(EEGEOGEMG)
% legend('EEG','EEG+EOG+EMG')


%% Unpaired t-test between human scoring and clustering (for EEG and EEGEOGEMG) (need to incorporate bonferroni correction)

load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');
load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/type1auc_human_scorers_table.mat');

% Remove LT duplicate
type1auc_human_scorers_table(86:90,:) = [];

idx_EEG = find(Table_AUC_100.NumChannels == 1);
Table_AUC_100 = Table_AUC_100(idx_EEG,:);

Subs = [1 5 439 458 596 748 749 752 604 807 821 870];

% Get mean across iterations for each dataset
for D = 1:12
    a = find(Table_AUC_100.Dataset == Subs(D));
    EEG(D) = mean(Table_AUC_100.AUC(a));    % AUC for EEG condition  
end

load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');
idx_EEGEOGEMG = find(Table_AUC_100.NumChannels == 3);
Table_AUC_100 = Table_AUC_100(idx_EEGEOGEMG,:);

for D = 1:12
    a = find(Table_AUC_100.Dataset == Subs(D));
    EEGEOGEMG(D) = mean(Table_AUC_100.AUC(a));  % AUC for EEG+EOG+EMG condition
end


%%%% Between novice and clustering, EEG condition
idx_nov_EEG = find(type1auc_human_scorers_table.channel == '1' & type1auc_human_scorers_table.dataset_type == 'Novice');
AUC_nov_EEG = str2double(type1auc_human_scorers_table.type1auc(idx_nov_EEG));

% Unpaired t-test
[h,p,ci,stats]= ttest2(AUC_nov_EEG,EEG)   % p = 0.0871   NOT SIGN

%%%% Between novice and clustering, EEG+EOG+EMG condition
idx_nov_EEGEOGEMG = find(type1auc_human_scorers_table.channel == '3' & type1auc_human_scorers_table.dataset_type == 'Novice');
AUC_nov_EEGEOGEMG = str2double(type1auc_human_scorers_table.type1auc(idx_nov_EEGEOGEMG));

% Unpaired t-test
[h,p,ci,stats]= ttest2(AUC_nov_EEGEOGEMG,EEGEOGEMG)   % p = 0.0137  clust SIGN higher

%%%% Between expert and clustering, EEG condition
idx_exp_EEG = find(type1auc_human_scorers_table.channel == '1' & type1auc_human_scorers_table.dataset_type == 'Expert');
AUC_exp_EEG = str2double(type1auc_human_scorers_table.type1auc(idx_exp_EEG));

% Unpaired t-test
[h,p,ci,stats]= ttest2(AUC_exp_EEG,EEG)   % p = 0.0516    NOT SIGN

%%%% Between expert and clustering, EEG+EOG+EMG condition
idx_exp_EEGEOGEMG = find(type1auc_human_scorers_table.channel == '3' & type1auc_human_scorers_table.dataset_type == 'Expert');
AUC_exp_EEGEOGEMG = str2double(type1auc_human_scorers_table.type1auc(idx_exp_EEGEOGEMG));

% Unpaired t-test
[h,p,ci,stats]= ttest2(AUC_exp_EEGEOGEMG,EEGEOGEMG)   % p = 0.0460   Cluster SIGN lower



%% Paired t-test; difference between the EEG and EOM for same group (not sign diff)

% Between EEG and EOM for expert scorers (need to remove LT duplicate!!)
[h,p,ci,stats] = ttest(AUC_exp_EEG,AUC_exp_EEGEOGEMG)  % X 

% Between EEG and EOM for novice scorers
[h,p,ci,stats] = ttest(AUC_nov_EEG,AUC_nov_EEGEOGEMG)  % X 

% Between EEG and EOM for clustering 
[h,p,ci,stats] = ttest(EEG,EEGEOGEMG)  % X 


%% Bonf test

X = [meanS1 meanP1; meanS2 meanP2; meanS3 meanP3; meanS4 meanP4; meanS5 meanP5];
[h,p,sigPairs] = ttest_bonf(X);

for D = 1:12
    L(D) = length(removed_feat_idx{1,D})
end

%% Table JASP

load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/type1auc_human_scorers_table.mat');
load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');

Subs = [1 5 439 458 596 748 749 752 604 807 821 870];

% Remove LT duplicate
type1auc_human_scorers_table(86:90,:) = [];

% AUC Cluster, EEG
idx_EEG = find(Table_AUC_100.NumChannels == 1);
Table_AUC_100 = Table_AUC_100(idx_EEG,:);

for D = 1:12
    a = find(Table_AUC_100.Dataset == Subs(D));
    AUC_clust_EEG(D) = mean(Table_AUC_100.AUC(a));      
end

% AUC Cluster, EOM
load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');
idx_EEGEOGEMG = find(Table_AUC_100.NumChannels == 3);
Table_AUC_100 = Table_AUC_100(idx_EEGEOGEMG,:);

for D = 1:12
    a = find(Table_AUC_100.Dataset == Subs(D));
    AUC_clust_EOM(D) = mean(Table_AUC_100.AUC(a));  
end

% AUC Novice, EEG 
idx_nov_EEG = find(type1auc_human_scorers_table.channel == '1' & type1auc_human_scorers_table.dataset_type == 'Novice');
AUC_nov_EEG = str2double(type1auc_human_scorers_table.type1auc(idx_nov_EEG));

    % Average the 5 stages for each subject
    S  = numel(AUC_nov_EEG);
    xx = reshape(AUC_nov_EEG(1:S - mod(S, 5)), 5, []);
    AUC_nov_EEG  = sum(xx, 1).' / 5;

% AUC Novice, EOM 
idx_nov_EEGEOGEMG = find(type1auc_human_scorers_table.channel == '3' & type1auc_human_scorers_table.dataset_type == 'Novice');
AUC_nov_EOM = str2double(type1auc_human_scorers_table.type1auc(idx_nov_EEGEOGEMG));

    % Average the 5 stages for each subject
    S  = numel(AUC_nov_EOM);
    xx = reshape(AUC_nov_EOM(1:S - mod(S, 5)), 5, []);
    AUC_nov_EOM  = sum(xx, 1).' / 5;

% AUC Expert, EEG 
idx_exp_EEG = find(type1auc_human_scorers_table.channel == '1' & type1auc_human_scorers_table.dataset_type == 'Expert');
AUC_exp_EEG = str2double(type1auc_human_scorers_table.type1auc(idx_exp_EEG));

    % Average the 5 stages for each subject
    S  = numel(AUC_exp_EEG);
    xx = reshape(AUC_exp_EEG(1:S - mod(S, 5)), 5, []);
    AUC_exp_EEG  = sum(xx, 1).' / 5;
    
% AUC Expert, EOM 
idx_exp_EEGEOGEMG = find(type1auc_human_scorers_table.channel == '3' & type1auc_human_scorers_table.dataset_type == 'Expert');
AUC_exp_EOM = str2double(type1auc_human_scorers_table.type1auc(idx_exp_EEGEOGEMG));

    % Average the 5 stages for each subject
    S  = numel(AUC_exp_EOM);
    xx = reshape(AUC_exp_EOM(1:S - mod(S, 5)), 5, []);
    AUC_exp_EOM  = sum(xx, 1).' / 5;

% Column with AUC EEG from Cluster, Novice and Expert
AUC_EEG_EOM = [AUC_clust_EEG' AUC_clust_EOM';AUC_nov_EEG AUC_nov_EOM;AUC_exp_EEG AUC_exp_EOM];
Table1 = array2table(AUC_EEG_EOM,'VariableNames',{'EEG','EOM'});

Dataset = {'1','2','3','4','5','6','7','8','9','10','11','12','1|2','1|2','1|2','1|2','1|2','1|2','1|2','1|2','1|2','1|2','1|2','1|2','1|2'}';
ScoringEntity = {'Clust','Clust','Clust','Clust','Clust','Clust','Clust','Clust','Clust','Clust','Clust','Clust','NOV','NOV','NOV','NOV','NOV','NOV','NOV','NOV','EXP','EXP','EXP','EXP','EXP'}';
Table2 = array2table([Dataset ScoringEntity],'VariableNames',{'Dataset','ScoringEntity'});

% Table
Table = [Table2 Table1];
fpath = '/Users/nico/Documents/HCTSA/Analysis/AUC_100';
% writetable(Table, [fpath filesep 'Table_AUC.csv'])  % save

%  Unpaired t test (difference EEG EOM)
[h,p,ci,stats] = ttest(Table.EEG(1:12,:),Table.EOM(1:12,:))  % X 


