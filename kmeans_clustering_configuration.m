output_folder=datestr(now,30);

addpath('library/distinguishable_colors');

% KJ_N1
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/KJ_N1/30seconds/HCTSA_N.mat';
% C4_COL=1:789;
% F4_COL=789*2+1:789*3;
% C4_O2_COL=1:789*2;
% C4_F4_COL=1:789*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=6;
% timeseries_sampling_rate=200;

% KJ_N1_bipolar
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/KJ_N1_bipolar/HCTSA_N.mat';
% C4_COL=1:789;
% EOG_COL=(length(C4_COL)+1):((length(C4_COL)+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*2+1):((length(C4_COL)*2+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=3;
% timeseries_sampling_rate=200;


% % KJ_N2
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/KJ_N2/HCTSA_N.mat';
% C4_COL=1:1322;
% C4_O2_COL=1:length(C4_COL)*2;
% C4_F4_COL=1:length(C4_COL)*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=6;
% timeseries_sampling_rate=200;

% % KJ_N2_bipolar
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/KJ_N2_bipolar/HCTSA_N.mat';
% C4_COL=1:1322;
% EOG_COL=(length(C4_COL)+1):((length(C4_COL)+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*2+1):((length(C4_COL)*2+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=3;
% timeseries_sampling_rate=200;

% LI_N2
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/LI_N2/HCTSA_N.mat';
% C4_COL=1:1375;
% C4_O2_COL=1:length(C4_COL)*2;
% C4_F4_COL=1:length(C4_COL)*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=6;
% timeseries_sampling_rate=200;

% LI_N2_bipolar
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/LI_N2_bipolar/HCTSA_N.mat';
% C4_COL=1:1375;
% EOG_COL=(length(C4_COL)+1):((length(C4_COL)+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*2+1):((length(C4_COL)*2+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=3;
% timeseries_sampling_rate=200;

% LI_N2_C4_C3
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/LI_N2_C4_C3/HCTSA_N.mat';
% C4_COL=1:1375;
% C4_O2_COL=1:length(C4_COL)*2;
% C4_F4_COL=1:length(C4_COL)*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=1;
% timeseries_sampling_rate=200;

% % ME_N1
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/ME_N1/HCTSA_N.mat';
% C4_COL=1:1120;
% C4_O2_COL=1:length(C4_COL)*2;
% C4_F4_COL=1:length(C4_COL)*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=6;
% timeseries_sampling_rate=200;

% % ME_N1_bipolar
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/ME_N1_bipolar/HCTSA_N.mat';
% C4_COL=1:1120;
% EOG_COL=(length(C4_COL)+1):((length(C4_COL)+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*2+1):((length(C4_COL)*2+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=3;
% timeseries_sampling_rate=200;

% % ME_N1_C4_C3
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/ME_N1_C4_C3/HCTSA_N.mat';
% C4_COL=1:1120;
% C4_O2_COL=1:length(C4_COL)*2;
% C4_F4_COL=1:length(C4_COL)*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=1;
% timeseries_sampling_rate=200;

% ME_N2
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/ME_N2/HCTSA_N.mat';
% C4_COL=1:1192;
% C4_O2_COL=1:length(C4_COL)*2;
% C4_F4_COL=1:length(C4_COL)*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=6;
% timeseries_sampling_rate=200;

% % ME_N2_bipolar
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/ME_N2_bipolar/HCTSA_N.mat';
% C4_COL=1:1192;
% EOG_COL=(length(C4_COL)+1):((length(C4_COL)+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*2+1):((length(C4_COL)*2+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=3;
% timeseries_sampling_rate=200;

% % ME_N3
% hctsafile = '/Volumes/Spaceship/Voss_Lucid/ME_N3/HCTSA_N.mat';
% C4_COL=1:957;
% C4_O2_COL=1:length(C4_COL)*2;
% C4_F4_COL=1:length(C4_COL)*3;
% EOG_COL=(length(C4_COL)*3+1):((length(C4_COL)*3+1)+length(C4_COL)-1);
% EMG_COL=(length(C4_COL)*5+1):((length(C4_COL)*5+1)+length(C4_COL)-1);
% epoch_seconds=30;
% no_of_channels=6;
% timeseries_sampling_rate=200;

% % ME_N3_bipolar
hctsafile = '/Volumes/Spaceship/Voss_Lucid/ME_N3_bipolar/HCTSA_N.mat';
C4_COL=1:957;
EOG_COL=(length(C4_COL)+1):((length(C4_COL)+1)+length(C4_COL)-1);
EMG_COL=(length(C4_COL)*2+1):((length(C4_COL)*2+1)+length(C4_COL)-1);
epoch_seconds=30;
no_of_channels=3;
timeseries_sampling_rate=200;

