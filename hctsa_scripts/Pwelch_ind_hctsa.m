

%% SET PARAMETERS

format long g

% savedir = '/Users/nico/Documents/Thesis Internship/WISC/Work WISC/Decat LaBerge Baird project/Material/LD_REMwoArous/';
rootdir = '/Users/nico/Documents/Thesis Internship/WISC/Work WISC/Decat LaBerge Baird project/Material/LD_REMwoArous/';

subs = {'eh0527dl','js1207el2c','nl1010gl','nl0530el','nl1013bl','nl1123b','eh0708jl','nl0530gl'};                                                                    % WILDs N2

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

% 
% if equalize_remlucidsegs == 1
%     if size(rem_epoch,2) > size(lucid_epoch,2)
%         rem_epoch = rem_epoch(:,end-size(lucid_epoch,2)+1:end);
%     elseif size(lucid_epoch,2) > size(rem_epoch,2)
%         lucid_epoch = lucid_epoch(:,1:size(rem_epoch,2));
%     end
% end


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
fft_lucid_deltaCTL(:,iFile) = squeeze(mean(fft_rem(:,delta_idx),2));
fft_lucid_thetaCTL(:,iFile) = squeeze(mean(fft_rem(:,theta_idx),2));
fft_lucid_alphaCTL(:,iFile) = squeeze(mean(fft_rem(:,alpha_idx),2));
fft_lucid_betaCTL(:,iFile) = squeeze(mean(fft_rem(:,beta_idx),2));
fft_lucid_gammaCTL(:,iFile) = squeeze(mean(fft_rem(:,gamma_idx),2));

% fft_rem_allfreqs = [fft_rem_delta fft_rem_theta fft_rem_alpha fft_rem_beta fft_rem_gamma];
% fft_rem_allsubs = cat(3,fft_rem_allsubs,fft_rem_allfreqs);

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
fft_lucid_delta(:,iFile) = squeeze(mean(fft_lucid(:,delta_idx),2));
fft_lucid_theta(:,iFile) = squeeze(mean(fft_lucid(:,theta_idx),2));
fft_lucid_alpha(:,iFile) = squeeze(mean(fft_lucid(:,alpha_idx),2));
fft_lucid_beta(:,iFile) = squeeze(mean(fft_lucid(:,beta_idx),2));
fft_lucid_gamma(:,iFile) = squeeze(mean(fft_lucid(:,gamma_idx),2));

% fft_lucid_allfreqs = [fft_lucid_delta fft_lucid_theta fft_lucid_alpha fft_lucid_beta fft_lucid_gamma];
% fft_lucid_allsubs = cat(3,fft_lucid_allsubs,fft_lucid_allfreqs);


%% PLOT

% generate individual topoplots


% saveas(gcf,sprintf('%s/%s_NonLucidSeg(30).fig',savedir,sub))
% saveas(gcf,sprintf('%s/%s_NonLucidSeg(30).jpg',savedir,sub))
% 
% close(gcf)


end

% Average power at each electrode across all lucid segments

for i = 1:28
    MeanSubsDelta(i,1) = mean(fft_lucid_delta(i,:));
    MeanSubsTheta(i,1) = mean(fft_lucid_theta(i,:));
    MeanSubsAlpha(i,1) = mean(fft_lucid_alpha(i,:));
    MeanSubsBeta(i,1) = mean(fft_lucid_beta(i,:));
    MeanSubsGamma(i,1) = mean(fft_lucid_gamma(i,:));
end

for i = 1:28
    MeanSubsDeltaCTL(i,1) = mean(fft_lucid_deltaCTL(i,:));
    MeanSubsThetaCTL(i,1) = mean(fft_lucid_thetaCTL(i,:));
    MeanSubsAlphaCTL(i,1) = mean(fft_lucid_alphaCTL(i,:));
    MeanSubsBetaCTL(i,1) = mean(fft_lucid_betaCTL(i,:));
    MeanSubsGammaCTL(i,1) = mean(fft_lucid_gammaCTL(i,:));
end



figure;
suptitle('Lucid segment average power (LR2 to awakening)')
subplot(1,5,1); topoplot(MeanSubsDelta, EEG.chanlocs(1:28),'maplimits',[0 10],'electrodes','on','style','both','headrad','rim'); title('Delta (2-4 Hz)'); colorbar;
subplot(1,5,2); topoplot(MeanSubsTheta, EEG.chanlocs(1:28),'maplimits',[0 5],'electrodes','on','style','both','headrad','rim'); title('Theta (5-8 Hz)'); colorbar;
subplot(1,5,3); topoplot(MeanSubsAlpha, EEG.chanlocs(1:28),'maplimits',[0 5],'electrodes','on','style','both','headrad','rim'); title('Alpha (8-12 Hz)'); colorbar;
subplot(1,5,4); topoplot(MeanSubsBeta, EEG.chanlocs(1:28),'maplimits',[0 0.3],'electrodes','on','style','both','headrad','rim'); title('Beta (15-35 Hz)'); colorbar;
subplot(1,5,5); topoplot(MeanSubsGamma, EEG.chanlocs(1:28),'maplimits',[0 0.08],'electrodes','on','style','both','headrad','rim'); title('Gamma (36-45 Hz)'); colorbar;
set(gcf,'Position',[0 0 1500 500])
set(gcf,'color','w');


figure;
suptitle('Nonlucid segment average power (LR2 to awakening)')
subplot(1,5,1); topoplot(MeanSubsDeltaCTL, EEG.chanlocs(1:28),'maplimits',[0 10],'electrodes','on','style','both','headrad','rim'); title('Delta (2-4 Hz)'); colorbar;
subplot(1,5,2); topoplot(MeanSubsThetaCTL, EEG.chanlocs(1:28),'maplimits',[0 5],'electrodes','on','style','both','headrad','rim'); title('Theta (5-8 Hz)'); colorbar;
subplot(1,5,3); topoplot(MeanSubsAlphaCTL, EEG.chanlocs(1:28),'maplimits',[0 5],'electrodes','on','style','both','headrad','rim'); title('Alpha (8-12 Hz)'); colorbar;
subplot(1,5,4); topoplot(MeanSubsBetaCTL, EEG.chanlocs(1:28),'maplimits',[0 0.3],'electrodes','on','style','both','headrad','rim'); title('Beta (15-35 Hz)'); colorbar;
subplot(1,5,5); topoplot(MeanSubsGammaCTL, EEG.chanlocs(1:28),'maplimits',[0 0.08],'electrodes','on','style','both','headrad','rim'); title('Gamma (36-45 Hz)'); colorbar;
set(gcf,'Position',[0 0 1500 500])
set(gcf,'color','w');


% saveas(gcf,'avg_N2_NonLucidSeg(30).fig');
% saveas(gcf,'avg_N2_NonLucidSeg(30).jpg');




