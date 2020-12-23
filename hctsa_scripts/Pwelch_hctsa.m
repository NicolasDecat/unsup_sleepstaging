
%% SET PARAMETERS

format long g
equalize_remlucidsegs = 0; % if 1 then set rem base = length of lucid REM else 0

% savedir = '/Users/nico/Documents/Thesis Internship/WISC/Work WISC/Decat LaBerge Baird project/Material/LD_PhasicTonicwoArous/Try005';
rootdir = '/Users/nico/Documents/Thesis Internship/WISC/Work WISC/Decat LaBerge Baird project/Material/LD_REMwoArous/';
rootdirCTL = '/Users/nico/Documents/Thesis Internship/WISC/Work WISC/Decat LaBerge Baird project/Material/LD_PhasicTonicwoArous/WILDs_TonicREM';
condition = 'ICApruned 2Hz';

% if equalize_remlucidsegs == 0
%     savedir = '/Volumes/BairdHD/Projects/LucidDreaming/LucidDreamingStudies/Spike40hz2019/Analysis/Analysis_060420/ICA_Final/Results_UnequalSegs';
% elseif equalize_remlucidsegs == 1
%     savedir = '/Volumes/BairdHD/Projects/LucidDreaming/LucidDreamingStudies/Spike40hz2019/Analysis/Analysis_060420/ICA_Final/Results_EqualSegs';
% end

subs = {'eh0527dl','js1207el2c','nl1010gl','nl0530el','nl1013bl','nl1123b','eh0708jl','nl0530gl'};
    %,'js1207el2c','nl1010gl','nl0530el','nl1013bl','nl1123b','eh0708jl','nl0530gl'}; % ,'eh0527g','nl0523cl'};
% subs = {'eh0527g','nl0523cl'};                                                                   

fft_rem_allsubs = [];
fft_lucid_allsubs = [];

for iFile = 1:length(subs)
    
sub = subs{iFile};

filename_rem = sprintf('%s/%s_REMwoArous.set',rootdir,sub);
filenam_lucid = sprintf('%s/%s_LR2Aw.set',rootdir,sub);

% load REM file
[ALLEEG,EEG,CURRENTSET,ALLCOM] = eeglab;
EEG = pop_loadset('filename',filename_rem);
eeglab redraw

rem_epoch = EEG.data;

% load LUCID file
[ALLEEG,EEG,CURRENTSET,ALLCOM] = eeglab;
EEG = pop_loadset('filename',filenam_lucid);
eeglab redraw

lucid_epoch = EEG.data;

% pwelch parameters
sampfreq = EEG.srate;
mspersamp = round(1000 / sampfreq);
window = round(2000 / mspersamp);
overlap = round(1000 / mspersamp);
nfft = round(2000 / mspersamp);


if equalize_remlucidsegs == 1
    if size(rem_epoch,2) > size(lucid_epoch,2)
        rem_epoch = rem_epoch(:,end-size(lucid_epoch,2)+1:end);
    elseif size(lucid_epoch,2) > size(rem_epoch,2)
        lucid_epoch = lucid_epoch(:,1:size(rem_epoch,2));
    end
end


%% COMPUTE FFT POWER

% REM
fft_rem = [];
freq_rem = [];
for i_el = 1:size(rem_epoch,1)
    [fft_rem(i_el,:),freq_rem] = pwelch(rem_epoch(i_el,:),window,overlap,nfft,sampfreq);
end

% find freq ranges of interest
delta_idx = find(freq_rem >= 2 & freq_rem <= 4);
theta_idx = find(freq_rem >= 5 & freq_rem <= 8);
alpha_idx = find(freq_rem >= 8 & freq_rem <= 12);
beta_idx = find(freq_rem >= 15 & freq_rem <= 35);
gamma_idx = find(freq_rem >= 36 & freq_rem <= 45);

% avg FFT in freq ranges of interest
fft_rem_delta = squeeze(mean(fft_rem(:,delta_idx),2));
fft_rem_theta = squeeze(mean(fft_rem(:,theta_idx),2));
fft_rem_alpha = squeeze(mean(fft_rem(:,alpha_idx),2));
fft_rem_beta = squeeze(mean(fft_rem(:,beta_idx),2));
fft_rem_gamma = squeeze(mean(fft_rem(:,gamma_idx),2));

fft_rem_allfreqs = [fft_rem_delta fft_rem_theta fft_rem_alpha fft_rem_beta fft_rem_gamma];
fft_rem_allsubs = cat(3,fft_rem_allsubs,fft_rem_allfreqs);

% LD
fft_lucid = [];
freq_lucid = [];
for i_el = 1:size(lucid_epoch,1)
    [fft_lucid(i_el,:),freq_lucid] = pwelch(lucid_epoch(i_el,:),window,overlap,nfft,sampfreq);
end

% find freq ranges of interest
delta_idx = find(freq_lucid >= 2 & freq_lucid <= 4);
theta_idx = find(freq_lucid >= 5 & freq_lucid <= 8);
alpha_idx = find(freq_lucid >= 8 & freq_lucid <= 12);
beta_idx = find(freq_lucid >= 15 & freq_lucid <= 35);
gamma_idx = find(freq_lucid >= 36 & freq_lucid <= 45);

% avg FFT in freq ranges of interest
fft_lucid_delta = squeeze(mean(fft_lucid(:,delta_idx),2));
fft_lucid_theta = squeeze(mean(fft_lucid(:,theta_idx),2));
fft_lucid_alpha = squeeze(mean(fft_lucid(:,alpha_idx),2));
fft_lucid_beta = squeeze(mean(fft_lucid(:,beta_idx),2));
fft_lucid_gamma = squeeze(mean(fft_lucid(:,gamma_idx),2));


fft_lucid_allfreqs = [fft_lucid_delta fft_lucid_theta fft_lucid_alpha fft_lucid_beta fft_lucid_gamma];
fft_lucid_allsubs = cat(3,fft_lucid_allsubs,fft_lucid_allfreqs);


%% PLOT

% generate difference topoplots
fft_diff_delta = fft_lucid_allfreqs(:,1) - fft_rem_allfreqs(:,1);
fft_diff_theta = fft_lucid_allfreqs(:,2) - fft_rem_allfreqs(:,2);
fft_diff_alpha = fft_lucid_allfreqs(:,3) - fft_rem_allfreqs(:,3);
fft_diff_beta = fft_lucid_allfreqs(:,4) - fft_rem_allfreqs(:,4);
fft_diff_gamma = fft_lucid_allfreqs(:,5) - fft_rem_allfreqs(:,5);

figure;
suptitle(sprintf('%s Phasic vs Tonic REM b005 ',sub))
subplot(1,5,1); topoplot(fft_diff_delta, EEG.chanlocs(1:28),'maplimits',[-3 3],'electrodes','on','style','both','headrad','rim'); title('Delta (2-4 Hz)'); colorbar;
subplot(1,5,2); topoplot(fft_diff_theta, EEG.chanlocs(1:28),'maplimits',[-3 3],'electrodes','on','style','both','headrad','rim'); title('Theta (5-8 Hz)'); colorbar;
subplot(1,5,3); topoplot(fft_diff_alpha, EEG.chanlocs(1:28),'maplimits',[-3 3],'electrodes','on','style','both','headrad','rim'); title('Alpha (8-12 Hz)'); colorbar;
subplot(1,5,4); topoplot(fft_diff_beta, EEG.chanlocs(1:28),'maplimits',[-0.05 0.05],'electrodes','on','style','both','headrad','rim'); title('Beta (15-35 Hz)'); colorbar;
subplot(1,5,5); topoplot(fft_diff_gamma, EEG.chanlocs(1:28),'maplimits',[-0.01 0.01],'electrodes','on','style','both','headrad','rim'); title('Gamma (36-45 Hz)'); colorbar;
set(gcf,'Position',[0 0 1500 500])
set(gcf,'color','w');

% saveas(gcf,sprintf('%s/%s_PhasicvsTonicb005.fig',savedir,sub))
% saveas(gcf,sprintf('%s/%s_PhasicvsTonicb005.jpg',savedir,sub))
% 
close(gcf)


end


% Compute difference across conditions
fft_diff_delta_allfiles = fft_lucid_allsubs(:,1,:) - fft_rem_allsubs(:,1,:);
fft_diff_theta_allfiles = fft_lucid_allsubs(:,2,:) - fft_rem_allsubs(:,2,:);
fft_diff_alpha_allfiles = fft_lucid_allsubs(:,3,:) - fft_rem_allsubs(:,3,:);
fft_diff_beta_allfiles = fft_lucid_allsubs(:,4,:) - fft_rem_allsubs(:,4,:);
fft_diff_gamma_allfiles = fft_lucid_allsubs(:,5,:) - fft_rem_allsubs(:,5,:);

% Average of "difference across conditions" across groups
fft_diff_delta_allfiles = mean(fft_diff_delta_allfiles,3);
fft_diff_theta_allfiles = mean(fft_diff_theta_allfiles,3);
fft_diff_alpha_allfiles = mean(fft_diff_alpha_allfiles,3);
fft_diff_beta_allfiles = mean(fft_diff_beta_allfiles,3);
fft_diff_gamma_allfiles = mean(fft_diff_gamma_allfiles,3);


figure;
suptitle('avg Phasic vs Tonic REM b005');
subplot(1,5,1); topoplot(fft_diff_delta_allfiles, EEG.chanlocs(1:28),'maplimits',[-3 3],'electrodes','on','style','both','headrad','rim'); title('Delta (2-4 Hz)'); colorbar;
subplot(1,5,2); topoplot(fft_diff_theta_allfiles, EEG.chanlocs(1:28),'maplimits',[-3 3],'electrodes','on','style','both','headrad','rim'); title('Theta (5-8 Hz)'); colorbar;
subplot(1,5,3); topoplot(fft_diff_alpha_allfiles, EEG.chanlocs(1:28),'maplimits',[-3 3],'electrodes','on','style','both','headrad','rim'); title('Alpha (8-12 Hz)'); colorbar;
subplot(1,5,4); topoplot(fft_diff_beta_allfiles, EEG.chanlocs(1:28),'maplimits',[-0.05 0.05],'electrodes','on','style','both','headrad','rim'); title('Beta (15-35 Hz)'); colorbar;
subplot(1,5,5); topoplot(fft_diff_gamma_allfiles, EEG.chanlocs(1:28),'maplimits',[-0.01 0.01],'electrodes','on','style','both','headrad','rim'); title('Gamma (36-45 Hz)'); colorbar;
set(gcf,'Position',[0 0 1500 500])
set(gcf,'color','w');

% saveas(gcf, fullfile(savedir,'avg_WILDvsTonicb005'),'jpg');
% saveas(gcf, fullfile(savedir,'avg_WILDvsTonicb005'),'fig');







