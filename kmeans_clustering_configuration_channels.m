    %% 
    dataset_names=["KJ_N1", "KJ_N2", "LI_N2", "ME_N1", "ME_N2", "ME_N3"];
%     dataset=find(strcmp(dataset_names,"KJ_N1"));
%     data(6,:)=data(5,:)-data(6,:);

%     channels{find(strcmp(dataset_names,"KJ_N1"))}.EEG = [11, 12];
    %channels{find(strcmp(dataset_names,"KJ_N1"))}.EEG = [3, 4, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25];
    %data(11,:)=data(11,:)-data(12,:);
    channels{find(strcmp(dataset_names,"KJ_N1"))}.EEG = [4, 12, 16]; % F4, C4, O2
    channels{find(strcmp(dataset_names,"KJ_N1"))}.EOG = [9,10];
    channels{find(strcmp(dataset_names,"KJ_N1"))}.EMG = [6];

    % KJ_N2
%    channels{find(strcmp(dataset_names,"KJ_N2"))}.EEG = [9,10];
    channels{find(strcmp(dataset_names,"KJ_N2"))}.EEG = [6, 10, 14]; % F4, C4, O2
%    channels{find(strcmp(dataset_names,"KJ_N2"))}.EEG = [5,6,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23];
    channels{find(strcmp(dataset_names,"KJ_N2"))}.EOG = [7,8];
    channels{find(strcmp(dataset_names,"KJ_N2"))}.EMG = [3];

    % % LI_N2
%    channels{find(strcmp(dataset_names,"LI_N2"))}.EEG = [9,10];
    channels{find(strcmp(dataset_names,"LI_N2"))}.EEG = [8, 10, 14]; % F4, C4, O2
%    channels{find(strcmp(dataset_names,"LI_N2"))}.EEG = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23];
    channels{find(strcmp(dataset_names,"LI_N2"))}.EOG = [1,2];
    channels{find(strcmp(dataset_names,"LI_N2"))}.EMG = [3];

% % ME_N1
%    channels{find(strcmp(dataset_names,"ME_N1"))}.EEG = [9,10];
    channels{find(strcmp(dataset_names,"ME_N1"))}.EEG = [8, 10, 14]; % F4, C4, O2
%    channels{find(strcmp(dataset_names,"ME_N1"))}.EEG = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23];
    channels{find(strcmp(dataset_names,"ME_N1"))}.EOG = [1,2];
    channels{find(strcmp(dataset_names,"ME_N1"))}.EMG = [3];

    % % ME_N2
%    channels{find(strcmp(dataset_names,"ME_N2"))}.EEG = [7, 8]; 
    channels{find(strcmp(dataset_names,"ME_N2"))}.EEG = [6, 8, 12]; % F4, C4, O2
%    channels{find(strcmp(dataset_names,"ME_N2"))}.EEG = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21];
    channels{find(strcmp(dataset_names,"ME_N2"))}.EOG = [1,2];
    channels{find(strcmp(dataset_names,"ME_N2"))}.EMG = [3];
    
%    channels{find(strcmp(dataset_names,"ME_N3"))}.EEG = [9,10];
    channels{find(strcmp(dataset_names,"ME_N3"))}.EEG = [8, 10, 14]; %F4, C4, O2
%    channels{find(strcmp(dataset_names,"ME_N3"))}.EEG = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23];
    channels{find(strcmp(dataset_names,"ME_N3"))}.EOG = [26,27];
    channels{find(strcmp(dataset_names,"ME_N3"))}.EMG = [3];

    dataset=find(strcmp(dataset_names, SBJ_ID));
    
    chans=channels{dataset};
    
    data(chans.EMG, :) = data(chans.EMG,:) - data(chans.EMG+1, :);
    
