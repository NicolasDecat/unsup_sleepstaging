config.base_dir="/Volumes/Spaceship/Voss_Lucid/";
config.hctsa_dir="/Users/Zhao/SleepPsychoPhysics/Source/hctsa";
config.hctsa_reduced_ops_file='reduced_ops.txt';

% Experiment related
config.no_of_channels_used=1;
config.run_base_folder='/1_CHANS_TST';
config.sub_clusters_range=[2:7];
config.total_channels_in_hctsa_datamat=3;

% Subject related.
config.subject_ids=["KJ_N2"];
config.subject_secondary_id="_bipolar";

% Load hctsa
homedir = pwd;
cd(char(config.hctsa_dir))
startup
cd(homedir)
