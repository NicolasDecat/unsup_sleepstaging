
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

% Average and stf for both conditions
mean_EEG = mean(EEG);
std_EEG = std(EEG);

mean_EEGEOGEMG = mean(EEGEOGEMG);
std_EEGEOGEMG = std(EEGEOGEMG);

% Paired t-test
[h,p,ci,stats] = ttest(EEGEOGEMG,EEG)

figure; plot(EEG)
hold on 
plot(EEGEOGEMG)
legend('EEG','EEG+EOG+EMG')

%% Unpaired t-test between human scoring and clustering (for EEG and EEGEOGEMG) (need to incorporate bonferroni correction)

load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/Table_AUC_100.mat');
load('/Users/nico/Documents/HCTSA/Analysis/AUC_100/type1auc_human_scorers_table.mat');

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
[h,p,ci,stats]= ttest2(AUC_exp_EEGEOGEMG,EEGEOGEMG)   % p = 0.0460   Cluster SIGN higher



%% In the paper, he got the average mean for each subject, then do paired t test again

% for example, try experts, EEG only
S1 = idx_exp_EEG(1:5);
S2 = idx_exp_EEG(6:10);
S3 = idx_exp_EEG(11:15);
S4 = idx_exp_EEG(16:20);
S5 = idx_exp_EEG(21:25);
S6 = idx_exp_EEG(26:30);

meanS1 = mean(str2double(type1auc_human_scorers_table.type1auc(S1)));
meanS2 = mean(str2double(type1auc_human_scorers_table.type1auc(S2)));
meanS3 = mean(str2double(type1auc_human_scorers_table.type1auc(S3)));
meanS4 = mean(str2double(type1auc_human_scorers_table.type1auc(S4)));
meanS5 = mean(str2double(type1auc_human_scorers_table.type1auc(S5)));
meanS6 = mean(str2double(type1auc_human_scorers_table.type1auc(S6)));

mean([meanS1 meanS2 meanS3 meanS4 meanS5 meanS6])
std([meanS1 meanS2 meanS3 meanS4 meanS5 meanS6])

% try experts, EEG EOG EMG

P1 = idx_exp_EEGEOGEMG(1:5);
P2 = idx_exp_EEGEOGEMG(6:10);
P3 = idx_exp_EEGEOGEMG(11:15);
P4 = idx_exp_EEGEOGEMG(16:20);
P5 = idx_exp_EEGEOGEMG(21:25);

meanP1 = mean(str2double(type1auc_human_scorers_table.type1auc(P1)));
meanP2 = mean(str2double(type1auc_human_scorers_table.type1auc(P2)));
meanP3 = mean(str2double(type1auc_human_scorers_table.type1auc(P3)));
meanP4 = mean(str2double(type1auc_human_scorers_table.type1auc(P4)));
meanP5 = mean(str2double(type1auc_human_scorers_table.type1auc(P5)));

mean([meanP1 meanP2 meanP3 meanP4 meanP5])
std([meanP1 meanP2 meanP3 meanP4 meanP5])


%%%% Bonf test

X = [meanS1 meanP1; meanS2 meanP2; meanS3 meanP3; meanS4 meanP4; meanS5 meanP5];

[h,p,sigPairs] = ttest_bonf(X);


for D = 1:12
    L(D) = length(removed_feat_idx{1,D})
end

mean(L)