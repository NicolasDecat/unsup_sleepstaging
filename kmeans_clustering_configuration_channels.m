    %% 
    dataset_names=["KJ_N1", "KJ_N2", "LI_N2", "ME_N1", "ME_N2"];
%     dataset=find(strcmp(dataset_names,"KJ_N1"));
%     data(6,:)=data(5,:)-data(6,:);

    channels{find(strcmp(dataset_names,"KJ_N1"))}.EEG = [11, 12];
    channels{find(strcmp(dataset_names,"KJ_N1"))}.EOG = [9,10];
    channels{find(strcmp(dataset_names,"KJ_N1"))}.EMG = [6];

    % KJ_N2
    % dataset=find(strcmp(dataset_names,"KJ_N2"));
    % data(9,:)=data(9,:)-data(24,:);
    % data(10,:)=data(10,:)-data(25,:);
    % data(3,:)=data(3,:)-data(4,:);
    channels{find(strcmp(dataset_names,"KJ_N2"))}.EEG = [9,10];
    channels{find(strcmp(dataset_names,"KJ_N2"))}.EOG = [7,8];
    channels{find(strcmp(dataset_names,"KJ_N2"))}.EMG = [3];

    % % LI_N2
%     dataset=find(strcmp(dataset_names,"LI_N2"));
%     data(3,:)=data(3,:)-data(4,:);
    channels{find(strcmp(dataset_names,"LI_N2"))}.EEG = [9,10];
    channels{find(strcmp(dataset_names,"LI_N2"))}.EOG = [1,2];
    channels{find(strcmp(dataset_names,"LI_N2"))}.EMG = [3];

% % ME_N1
%     dataset=find(strcmp(dataset_names,"ME_N1"));
%     data(3,:)=data(3,:)-data(4,:);
    channels{find(strcmp(dataset_names,"ME_N1"))}.EEG = [9,10];
%     channels{find(strcmp(dataset_names,"ME_N1"))}.EEG = [23];
    channels{find(strcmp(dataset_names,"ME_N1"))}.EOG = [1,2];
    channels{find(strcmp(dataset_names,"ME_N1"))}.EMG = [3];

    % % ME_N2
    dataset=find(strcmp(dataset_names,"ME_N2"));
    data(3,:)=data(3,:)-data(4,:);
    channels{find(strcmp(dataset_names,"ME_N2"))}.EEG = [7,8];
    channels{find(strcmp(dataset_names,"ME_N2"))}.EOG = [1,2];
    channels{find(strcmp(dataset_names,"ME_N2"))}.EMG = [3];
    
    chans=channels{dataset};
    