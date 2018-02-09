% xmldoc = xmlread(fullfile('learn-nsrr01-profusion.xml'))
% 
% xmlwrite(xmldoc)

addpath(genpath('C:\Users\Piengkwan\Documents\MATLAB\unsup_sleep_staging'))

data = xml2struct('learn-nsrr03-profusion.xml') % From http://au.mathworks.com/matlabcentral/fileexchange/28518-xml2struct

%% SleepStages 
% Epoch length
epochLength = data.CMPStudyConfig.EpochLength;

% Sleep stage data
sleepstageS = data.CMPStudyConfig.SleepStages.SleepStage;

%% Extract sleep stage in number
for i=1:length(sleepstageS)
    sleepstage(i,1) = str2num(sleepstageS{1,i}.Text);
end

%% 
plot(sleepstage)
axis([0 length(sleepstage) 0 6])

%% Save
save('Learn03Annot','sleepstage','epochLength')