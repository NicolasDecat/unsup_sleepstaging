%% Configuration
addpath('/Users/Zhao/Documents/MATLAB/Add-Ons/Collections/EDAToolboxV2');

% SOURCE_FILE='/Volumes/Spaceship/Voss_Lucid/ME_N1/ALL_EEG/HCTSA_N_ME_N1_1_EEG_Main.mat';
% TARGET_FILE='/Volumes/Spaceship/Voss_Lucid/ME_N2/ALL_EEG/HCTSA_N_ME_N2_1_EEG_Main.mat';

SOURCE_FILE='/Volumes/Spaceship/Voss_Lucid/ME_N1/ALL_EEG/HCTSA_N_ME_N1_1_EEG_7_REM_substages.mat';
TARGET_FILE='/Volumes/Spaceship/Voss_Lucid/ME_N2/ALL_EEG/HCTSA_N_ME_N2_1_EEG_7_REM_substages.mat';

homedir=pwd;

%%
t1=load(SOURCE_FILE);
ts1=struct2table(t1.TimeSeries);
to1=t1.Operations;
datamat1=t1.TS_DataMat;
ts1_labels = table2array(ts1(:,1))';
ts1_algo_groups = str2num(char(table2array(ts1(:,2))));
unique_groups1=unique(ts1_algo_groups);  

t2=load(TARGET_FILE);
ts2=struct2table(t2.TimeSeries);
to2=t2.Operations;
datamat2=t2.TS_DataMat;
ts2_labels = table2array(ts2(:,1))';
ts2_algo_groups = str2num(char(table2array(ts2(:,2))));
unique_groups2=unique(ts2_algo_groups);  

%% Synchronise the features (get the intersection)

OPS_FILE='reduced_ops.txt';
%%
fileID = fopen(OPS_FILE);
features = textscan(fileID,'%s %s %s');
fclose(fileID);

%% Wanted operation names
feat_name = features{1,2};

%% Check operation name, get feat_id
nn=0;
for n = 1:length(feat_name)
    op_name = char(feat_name(n));
    nn=nn+1;
    feat_id(nn) = nn;
    feat(nn).id = nn; % all_op.Operations(i).ID % Actual operation ID
    feat(nn).name = op_name;
end
clear i n nn op_name name

%% Loop through the top 200 operations and find out which features are missing in one but not the other
feat1_toremove = [];
feat2_toremove = [];

for i = 1:length(feat_name)
   curr_feat = char(feat_name(i));
   
   feat1_index = 0;
   feat2_index = 0;
   
   for f1 = 1:size(to1, 1)
      if strcmp(to1(f1).Name, curr_feat)
          feat1_index = f1;
          break;
      end
   end
   
   for f2 = 1:size(to2, 1)
      if strcmp(to2(f2).Name, curr_feat)
          feat2_index = f2;
          break;
      end
   end
   
   if (feat1_index == 0 && feat2_index ==0)
       disp(strcat('Both dataset 1 and 2 have no feature named ', curr_feat));
   elseif (feat1_index == 0)
       disp(strcat('Dataset 1 has no feature named: ', curr_feat, '. Going to remove from Dataset 2'));
       feat2_toremove = [feat2_toremove feat2_index];
   elseif (feat2_index == 0)
       disp(strcat('Dataset 2 has no feature named: ', curr_feat, '. Going to remove from Dataset 1'));
       feat1_toremove = [feat1_toremove feat1_index];
   end
end

%%
feat1_indices=1:size(datamat1, 2);
feat2_indices=1:size(datamat2, 2);

if length(feat1_toremove) ~= 0
    feat1_indices(feat1_indices==feat1_toremove) = [];
end

if length(feat2_toremove) ~= 0
    feat2_indices(feat2_indices==feat2_toremove) = [];
end

datamat1=datamat1(:, feat1_indices);
datamat2=datamat2(:, feat2_indices);

%% Perform the correlation for each group

for i = 1:size(unique_groups1, 1)
    ts1_grp_index = find(ts1_algo_groups == i);
    ts1_grp = datamat1(ts1_grp_index, :);
    
    feature_corr=zeros(size(unique_groups2, 1), size(datamat1, 2));
    corrs=[];
    
    for j = 1:size(unique_groups2, 1)
        ts2_grp_index = find(ts2_algo_groups == j);
        ts2_grp = datamat2(ts2_grp_index, :);
        
        t1g = ts1_grp;
        t2g = ts2_grp;
        
        if size(t1g, 1) > size(t2g, 1)
            rnd_indexes = randperm(size(t2g, 1));
            t1g = t1g(rnd_indexes, :);
        elseif size(t1g, 1) < size(t2g, 1)
            rnd_indexes = randperm(size(t1g, 1));
            t2g = t2g(rnd_indexes, :);
        end
       
        disp(strcat('Calculating adjust ran for ', num2str(i), '-', num2str(j)));
        R = corrcoef(t1g, t2g);
        corrs = [corrs R(1,2)];
        % Perform correlation for each feature (assume sorted)
%         for k = 1:size(datamat1, 2)
%             dm1 = ts1_grp(:,k);
%             dm2 = ts2_grp(:,k);
%             
%             % Check length and make them equal for correlation
%             if size(dm1, 1) > size(dm2, 1)
%                 rnd_indexes = randperm(size(dm2, 1));
%                 dm1 = dm1(rnd_indexes, :);
%             elseif size(dm1, 1) < size(dm2, 1)
%                 rnd_indexes = randperm(size(dm1, 1));
%                 dm2 = dm2(rnd_indexes, :);
%             end
%             
%             R = corrcoef(dm1, dm2);
% %             centre_distance = sqrt(mean(dm1)^2 + mean(dm2)^2);
% %             [tbl, chi2stat, pval]=crosstab(dm1,dm2);
% %             adjrand(dm1, dm2)
% %             
% %             if (pval < 0.05)
% %                 feature_corr(j, k) = 0;
% %             else
% %                 feature_corr(j, k) = 1;
% %             end
%             
%             feature_corr(j, k) = R(1, 2);
%             disp(R(1,2));
%             disp(mean(dm1));
%             disp(mean(dm2));
%                             
%             figure;
%             gscatter(dm1, dm2, '', 'rgb');
%             lsline;
%             xlim([0 1]);
%             ylim([0 1]);
%             close;
%         end
    end
    disp(strcat('------------ Dataset 1: Cluster ', num2str(i), '---------------'));
%     disp(mean(abs(feature_corr)'));
    disp(corrs);
end


