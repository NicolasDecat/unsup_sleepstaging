
%% Compare CCSHS, experts, Clust


%%%%% Let's focus on dataset 005 (3 experts in 001 but only 2 have same
%%%%% epochs)

%%%% CCSHS  (original labels)
load('/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/ccshs(3ch)/psychphysics_data')

Info = human_data(4271:end,1:2);
OrigLabels = human_data(4271:end,3);
Experts = human_data(4271:end,5);

%%%% Algorithm 

% 001: LT,SL,SN
% OO5: CD,PV   

% Get all cluster decisions _ Dataset 005
load('/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/algo(3ch)/statsOut_3ch_005')
ClustDecisions = statsOut.predictTest;
ClustDecisions_005 = changem(ClustDecisions,[1 2 3 4 5],[0 1 2 3 5]); % Convert [0 1 2 3 5] to [1 2 3 4 5]

%% Get cluster decision for each of the 5 experts, in order: CD,LT,PV,SL,SN (as in Info)
FileName = {'CD_1_3_5','PV_1_3_5'};
Experts_all = [];
EpochID = []; 
Conf = [];

%%%%% Let's focus on dataset 005 (3 experts in 001 but only 2 have same
%%%%% epochs)

% First, get scores from both experts, and from CCSHS
Expert_name = {"CD","PV"};
human_data.confidence = changem(human_data.confidence,[4 3 2 1],[1 2 3 4]);

for E = 1:2
    
    Expert_ID = find(human_data.subject == sprintf(Expert_name{E}) & (human_data.channels == '3'));
    Expert_label(:,E) = human_data.response(Expert_ID);
    CCSHS_label = human_data.label(Expert_ID);
    
    % To use later (confidence): flip confidence scale and get idx
    Conf = [Conf human_data.confidence(Expert_ID)];
    
end

Conf = mean(Conf')';  % Average confidence level across the 2 participants

% Then, get label from algorihm
for Subj = 1   % Whichever between CD an PV, cuz same epoch ID
    
    % Get the cluster decisions for one of the subjects
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/exp(3ch)/%s',FileName{Subj}),'trial')
    epoch_ID = [trial.fileID].';
    Algo_label = ClustDecisions_005(epoch_ID)';

end

Table = array2table([epoch_ID CCSHS_label Expert_label(:,1) Expert_label(:,2) Algo_label Conf],'VariableNames',{'EpochID','CCSHS','Expert1','Expert2','Algo','ExpertsConfidence'});

%%% 1) Instances where all entities agree: 51/205 (24.9%)

All_right = find(Table.CCSHS == Table.Expert1 & Table.Expert1 == Table.Expert2 & Table.Expert1 == Table.Algo);
Conf_All_right = Table.ExpertsConfidence(All_right);

%%% 2) Instances where humans agree but algo got something different 40/205 (19.5%)

Algo_weird = find(Table.CCSHS == Table.Expert1 & Table.Expert1 == Table.Expert2 & Table.Expert1 ~= Table.Algo);
Conf_Algo_weird = Table.ExpertsConfidence(Algo_weird);

%%% Disagreement within humans (55.6%) 
%%% 3a) Instances where algo got something different than  humans: 23/205 (11.2%)

Disa_diff=[];

for i = 1:205
     UNI = unique(table2array(Table(i,2:4)));
     if length(UNI) > 1 & Table.Algo(i) ~= UNI
        Disa_diff = [Disa_diff i];
     end
end

Conf_Disa_diff = Table.ExpertsConfidence(Disa_diff);

%%% 3b) Instances where algo got something like at least one humans: 91/205
%%% (44.4 - 3.9 = 40.1%)

Disa_same=[];

for i = 1:205
     UNI = unique(table2array(Table(i,2:4)));
     if length(UNI) > 1 & any(Table.Algo(i) == UNI)
        Disa_same = [Disa_same i];
     end
end

Conf_Disa_same = Table.ExpertsConfidence(Disa_same);

%%% 4) instances where ccshs got something different than experts and
%%% algo --> Effect of contextual information: 8/205 (3.9%) - included in
%%% the Case 3a)

CCSHS_diff = find(Table.Expert1 == Table.Expert2 & Table.Expert1 == Table.Algo & Table.Expert1 ~= Table.CCSHS);

Conf_CCSHS_diff = Table.ExpertsConfidence(CCSHS_diff);


%%% Table with all confidence: Conf_All_right, Conf_Algo_weird,
%%% Conf_disagr, Conf_CCSHS_diff

% Remove the 8 CCSHS_diff from Algo_Same
Disa_Same_noCCSHS = setdiff(Disa_same,CCSHS_diff);
Conf_Disa_Same_noCCSHS = Table.ExpertsConfidence(Disa_Same_noCCSHS);
Conf_disagr = [Conf_Disa_diff;Conf_Disa_Same_noCCSHS];

%%%% Conf_disagr SIGN different from Conf_all_right 
[h,p,ci,stats] = ttest2(Conf_All_right,Conf_disagr)


%%% Plot these instances: Algo_diff
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'

% Load TimeSeries
load('HCTSA_N.mat', 'TimeSeries')
LengthTS = size(TimeSeries,1)/7;
Raw_index = [1 2 3 4 5];  
Name_ST = {'W','N1','N2','N3','R'};    % Name of Cluster

Name = {'CW','C1','C2','C3','CR'};    % Name of Cluster
Equi_index = [0 1 2 3 5];             % Index of stages

%%% For given case, plot EEG+EOG+EMG

whichCase = CCSHS_diff;

% Get the epoch ID for corresponding case
NumTS = epoch_ID(whichCase);
    
% Plot
%%
for N = 1:length(NumTS)
    
    f = figure; ax=gca; 
    
    EEG = TimeSeries.Data(NumTS(N),:)+0.45;
    EOG = TimeSeries.Data(NumTS(N)+LengthTS,:)+0.25;
    EMG = TimeSeries.Data(NumTS(N)+LengthTS+LengthTS,:)+0.15;

    Channels = [{EEG} {EOG} {EMG}];
      
    for i = 1:3   % for each of the 3 channels (EEG, EOG and EMG - in one figure)

        plot(Channels{i},'LineWidth',1.5,'Color',[0, 0.4470, 0.7410])
        hold on
        
    end
    
    ax.Position = [0.09,0.23,0.89,0.65];
    
    % Xticks
    xticks([0:1280:3840]);
    ax.XTickLabels = {'0','','','30'};
    
    % Yticks
    yticks([0.15 0.30 0.45]);

    Labels = {'EMG','EOG','EEG'};
    ax.YAxis.TickLabels = Labels; % set
    % ax.YRuler.Axle.Visible = 'off'; 
    ax.YAxis.TickLength = [0 0];

    ylim([0 0.6])
    xlim([0 3840])

    % Get xlabel closter to axis
    xh = get(gca,'xlabel');
    p = get(xh,'position'); % get the current position property
    p(2) = -0.1 ;
    xlabel('Time (s)')
    set(xh,'position',p)   
    
    f.Position = [1,421,1440,376];
    set(gca,'box','off') 
    ax.FontSize = 25;
    
    Humans = strjoin([Name_ST(Table.CCSHS(whichCase(N))) Name_ST(Table.Expert1(whichCase(N))) Name_ST(Table.Expert2(whichCase(N)))]);
    Algo = Name(Table.Algo(whichCase(N)));
    
    title(sprintf('Humans: %s,  Algo: %s  (%s)',string(Humans), string(Algo), string(NumTS(N))))
    
end



%% Plot ind

for N = find(NumTS ==390)
    
    f = figure; ax=gca; 
    
    EEG = TimeSeries.Data(NumTS(N),:)+0.45;
    EOG = TimeSeries.Data(NumTS(N)+LengthTS,:)+0.25;
    EMG = TimeSeries.Data(NumTS(N)+LengthTS+LengthTS,:)+0.15;

    Channels = [{EEG} {EOG} {EMG}];
      
    for i = 1:3   % for each of the 3 channels (EEG, EOG and EMG - in one figure)

        plot(Channels{i},'LineWidth',1.5,'Color',[0, 0.4470, 0.7410])
        hold on
        
    end
    
    ax.Position = [0.09,0.23,0.89,0.65];
    
    % Xticks
    xticks([0:1280:3840]);
    ax.XTickLabels = {'0','','','30'};
    
    % Yticks
    yticks([0.15 0.30 0.45]);

    Labels = {'EMG','EOG','EEG'};
    ax.YAxis.TickLabels = Labels; % set
    %ax.YRuler.Axle.Visible = 'off'; 
    ax.YAxis.TickLength = [0 0];

    ylim([0 0.6])
    xlim([0 3840])

    % Get xlabel closter to axis
    xh = get(gca,'xlabel');
    p = get(xh,'position'); % get the current position property
    p(2) = -0.05 ;
    xlabel('Time (s)')
    set(xh,'position',p)   
    
    f.Position = [1,490,1392,307];
    set(gca,'box','off') 
    ax.FontSize = 25;
    
%     Humans = Name_ST(Table.Expert1(Algo_diff(N)));
%     Algo = Name(Table.Algo(Algo_diff(N)));
    
    %title(sprintf('Humans: %s,  Algo: %s  (%s)',string(Humans), string(Algo), string(NumTS(N))))
    
end

    
    fpath = '/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/Case3';
    set(f, 'Color', 'w')
    addpath '/Users/nico/Documents/MATLAB/hctsexport_fig-mastea-master/r'
   % export_fig([fpath filesep 'N1RRCR'],'-r 300')




%%% Try
% whichTS = 527;
% TS = TimeSeries.Data(whichTS,:);
% f = figure; plot(TS)
% f.Position = [1,608,1400,113];
% % For title
% title(sprintf('CCSHS: N2   Experts/Algo: REM',Stage{stage_human},Stage{stage_algo}))
%     

%%% Plot these instances: CCSHS_diff
