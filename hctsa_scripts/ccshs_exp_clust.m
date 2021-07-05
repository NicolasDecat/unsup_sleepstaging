
%% Compare CCSHS, experts, Clust  (005)

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
[h,p,ci,stats] = ttest2(Conf_All_right,Conf_disagr);


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

whichCase = Algo_weird;

% Get the epoch ID for corresponding case
NumTS = epoch_ID(whichCase);


%% Epochs i chose

%%% Case 1 (disagr)
% REM N2 REM CR: 732
% N3 N2 N3 C2: 461

%%% Case 2 (CCSHS diff)
% R N2 N2 C2: 1069
% N1 REM REM CR: 390

%%% Case 3 (Algo diff)
% W W W C1: 945
% N2 N2 N2 CR: 972
    
%% Plot

% for N = 1:length(NumTS)
%     
%     f = figure; ax=gca; 
%     
%     EEG = TimeSeries.Data(NumTS(N),:)+0.45;
%     EOG = TimeSeries.Data(NumTS(N)+LengthTS,:)+0.25;
%     EMG = TimeSeries.Data(NumTS(N)+LengthTS+LengthTS,:)+0.15;
% 
%     Channels = [{EEG} {EOG} {EMG}];
%       
%     for i = 1:3   % for each of the 3 channels (EEG, EOG and EMG - in one figure)
% 
%         plot(Channels{i},'LineWidth',1.5,'Color',[0, 0.4470, 0.7410])
%         hold on
%         
%     end
%     
%     ax.Position = [0.09,0.23,0.89,0.65];
%     
%     % Xticks
%     xticks([0:1280:3840]);
%     ax.XTickLabels = {'0','','','30'};
%     
%     % Yticks
%     yticks([0.15 0.30 0.45]);
% 
%     Labels = {'EMG','EOG','EEG'};
%     ax.YAxis.TickLabels = Labels; % set
%     % ax.YRuler.Axle.Visible = 'off'; 
%     ax.YAxis.TickLength = [0 0];
% 
%     ylim([0 0.6])
%     xlim([0 3840])
% 
%     % Get xlabel closter to axis
%     xh = get(gca,'xlabel');
%     p = get(xh,'position'); % get the current position property
%     p(2) = -0.1 ;
%     xlabel('Time (s)')
%     set(xh,'position',p)   
%     
%     f.Position = [1,421,1440,376];
%     set(gca,'box','off') 
%     ax.FontSize = 25;
%     
%     Humans = strjoin([Name_ST(Table.CCSHS(whichCase(N))) Name_ST(Table.Expert1(whichCase(N))) Name_ST(Table.Expert2(whichCase(N)))]);
%     Algo = Name(Table.Algo(whichCase(N)));
%     
%     title(sprintf('Humans: %s,  Algo: %s  (%s)',string(Humans), string(Algo), string(NumTS(N))))
%     
% end



%% Plot ind

for N = 972
    
    f = figure; ax=gca; 
    
    EEG = TimeSeries.Data((N),:)+0.45;
    EOG = TimeSeries.Data((N)+LengthTS,:)+0.25;
    EMG = TimeSeries.Data((N)+LengthTS+LengthTS,:)+0.15;

    Channels = [{EEG} {EOG} {EMG}];
      
    for i = 1:3   % for each of the 3 channels (EEG, EOG and EMG - in one figure)

        plot(Channels{i},'LineWidth',1.2,'Color',[0, 0.4470, 0.7410])
        hold on
        
    end
        
    % Xticks
    xticks([0:1280:3840]);
    ax.XTickLabels = {'0','','','30'};
    ax.XTickLabels = {'','','',''};
    
    % Yticks
    yticks([0.12 0.26 0.41]);

    Labels = {'EMG','EOG','EEG'};
    ax.YAxis.TickLabels = Labels; % set
    ax.YAxis.TickLength = [0 0];
 
    ylim([0 0.6])
    xlim([0 3840])

    % Get xlabel closter to axis
    xh = get(gca,'xlabel');
    p = get(xh,'position'); % get the current position property
    p(2) = -0.05 ;
    xlabel('')
    set(xh,'position',p)   
    
    f.Position = [1,490,725,307];
    ax.Position = [0.09,0.1,0.89,0.75];
        
    set(gca,'box','off') 
    ax.FontSize = 23;
    
%     Humans = Name_ST(Table.Expert1(Algo_diff(N)));
%     Algo = Name(Table.Algo(Algo_diff(N)));
    
    %title(sprintf('Humans: %s,  Algo: %s  (%s)',string(Humans), string(Algo), string(NumTS(N))))
     T = title('3ii');
     T.Position = [1920.002205782158,0.543156206054439,1.4210854715202e-14];
     
     ax.YRuler.Axle.Visible = 'off'; 

end


fpath = '/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/Case3';
set(f, 'Color', 'w')
addpath '/Users/nico/Documents/MATLAB/hctsexport_fig-mastea-master/r'
%export_fig([fpath filesep 'N2CR_sm'],'-r 300')


%%% Try
% whichTS = 527;
% TS = TimeSeries.Data(whichTS,:);
% f = figure; plot(TS)
% f.Position = [1,608,1400,113];
% % For title
% title(sprintf('CCSHS: N2   Experts/Algo: REM',Stage{stage_human},Stage{stage_algo}))
%     

%%% Plot these instances: CCSHS_diff






%% Now, focus on 001

%%%% CCSHS  (original labels)
load('/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/ccshs(3ch)/psychphysics_data')

Info = human_data(4271:end,1:2);
OrigLabels = human_data(4271:end,3);
Experts = human_data(4271:end,5);

%%%% Algorithm 

% 001: LT,SL,SN,   but remove LT cuz different epochs 

% Get all cluster decisions _ Dataset 001
load('/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/algo(3ch)/statsOut_3ch_001')
ClustDecisions = statsOut.predictTest;
ClustDecisions_001 = changem(ClustDecisions,[1 2 3 4 5],[0 1 2 3 5]); % Convert [0 1 2 3 5] to [1 2 3 4 5]

%%% Get cluster decision for each of the 5 experts, in order: CD,LT,PV,SL,SN (as in Info)
FileName = {'SL_2_3__','SN_1_3_1'};
Experts_all = [];
EpochID = []; 
Conf_001 = [];

%%%%% Let's focus on dataset 005 (3 experts in 001 but only 2 have same
%%%%% epochs)

% First, get scores from both experts, and from CCSHS
Expert_name = {"SL","SN"};
% human_data.confidence = changem(human_data.confidence,[4 3 2 1],[1 2 3
% 4]);   uncomment if Section 1 not ran

for E = 1:2
    
    Expert_ID_001 = find(human_data.subject == sprintf(Expert_name{E}) & (human_data.channels == '3'));
    Expert_label_001(:,E) = human_data.response(Expert_ID_001);
    CCSHS_label_001 = human_data.label(Expert_ID_001);
    
    % To use later (confidence): flip confidence scale and get idx
    Conf_001 = [Conf_001 human_data.confidence(Expert_ID_001)];
    
end

Conf_001 = mean(Conf_001')';  % Average confidence level across the 2 participants

% Then, get label from algorihm
for Subj = 1   % Whichever between SL an SN, cuz same epoch ID
    
    % Get the cluster decisions for one of the subjects
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo/exp(3ch)/%s',FileName{Subj}),'trial')
    epoch_ID_001 = [trial.fileID].';
    Algo_label_001 = ClustDecisions_001(epoch_ID_001)';

end

Table_001 = array2table([epoch_ID_001 CCSHS_label_001 Expert_label_001(:,1) Expert_label_001(:,2) Algo_label_001 Conf_001],'VariableNames',{'EpochID','CCSHS','Expert1','Expert2','Algo','ExpertsConfidence'});

%%% 1) Instances where all entities agree: 42/185 (22.7%)

All_right_001 = find(Table_001.CCSHS == Table_001.Expert1 & Table_001.Expert1 == Table_001.Expert2 & Table_001.Expert1 == Table_001.Algo);
Conf_All_right_001 = Table_001.ExpertsConfidence(All_right_001);

%%% 2) Instances where humans agree but algo got something different 60/185 (32.4%)

Algo_weird_001 = find(Table_001.CCSHS == Table_001.Expert1 & Table_001.Expert1 == Table_001.Expert2 & Table_001.Expert1 ~= Table_001.Algo);
Conf_Algo_weird_001 = Table_001.ExpertsConfidence(Algo_weird_001);

%%% Disagreement within humans (44.3%) 
%%% 3a) Instances where algo got something different than  humans: 30/185 (16.2%)

Disa_diff_001=[];

for i = 1:185
     UNI = unique(table2array(Table_001(i,2:4)));
     if length(UNI) > 1 & Table_001.Algo(i) ~= UNI
        Disa_diff_001 = [Disa_diff_001 i];
     end
end

Conf_Disa_diff_001 = Table_001.ExpertsConfidence(Disa_diff_001);

%%% 3b) Instances where algo got something like at least one humans: 53/185
%%% (28.6 - 0.6 = 28%)

Disa_same_001=[];

for i = 1:185
     UNI = unique(table2array(Table_001(i,2:4)));
     if length(UNI) > 1 & any(Table_001.Algo(i) == UNI)
        Disa_same_001 = [Disa_same_001 i];
     end
end

Conf_Disa_same_001 = Table_001.ExpertsConfidence(Disa_same_001);

%%% 4) instances where ccshs got something different than experts and
%%% algo --> Effect of contextual information: 1/185 (0.6%) - included in
%%% the Case 3a)

CCSHS_diff_001 = find(Table_001.Expert1 == Table_001.Expert2 & Table_001.Expert1 == Table_001.Algo & Table_001.Expert1 ~= Table_001.CCSHS);

Conf_CCSHS_diff_001 = Table_001.ExpertsConfidence(CCSHS_diff_001);


%%% Table with all confidence: Conf_All_right, Conf_Algo_weird,
%%% Conf_disagr, Conf_CCSHS_diff

% Remove the 8 CCSHS_diff from Algo_Same
Disa_Same_noCCSHS_001 = setdiff(Disa_same_001,CCSHS_diff_001);
Conf_Disa_Same_noCCSHS_001 = Table_001.ExpertsConfidence(Disa_Same_noCCSHS_001);
Conf_disagr_001 = [Conf_Disa_diff_001;Conf_Disa_Same_noCCSHS_001];


%% Compute agreement between all entities

% Exp1 vs. Exp2
A_expexp = find(Expert_label(:,1) == Expert_label(:,2));
A_expexp_perc = numel(A_expexp)/205*100;

A_expexp_001 = find(Expert_label_001(:,1) == Expert_label_001(:,2));
A_expexp_001_perc = numel(A_expexp_001)/185*100;

A_expexp_PERC = mean([A_expexp_perc A_expexp_001_perc]);

% Exp1 vs AASM
A_exp1AASM = find(Expert_label(:,1) == CCSHS_label);
A_exp1AASM_perc = numel(A_exp1AASM)/205*100;

A_exp1AASM_001 = find(Expert_label_001(:,1) == CCSHS_label_001);
A_exp1AASM_001_perc = numel(A_exp1AASM_001)/185*100;

A_exp1AASM_PERC = mean([A_exp1AASM_perc A_exp1AASM_001_perc]);

% Exp2 vs AASM
A_exp2AASM = find(Expert_label(:,2) == CCSHS_label);
A_exp2AASM_perc = numel(A_exp2AASM)/205*100;

A_exp2AASM_001 = find(Expert_label_001(:,2) == CCSHS_label_001);
A_exp2AASM_001_perc = numel(A_exp2AASM_001)/185*100;

A_exp2AASM_PERC = mean([A_exp2AASM_perc A_exp2AASM_001_perc]);

% Exp1 vs Algo
A_exp1Algo = find(Expert_label(:,1) == Algo_label);
A_exp1Algo_perc = numel(A_exp1Algo)/205*100;

A_exp1Algo_001 = find(Expert_label_001(:,1) == Algo_label_001);
A_exp1Algo_001_perc = numel(A_exp1Algo_001)/185*100;

A_exp1Algo_PERC = mean([A_exp1Algo_perc A_exp1Algo_001_perc]);

% Exp2 vs Algo
A_exp2Algo = find(Expert_label(:,2) == Algo_label);
A_exp2Algo_perc = numel(A_exp2Algo)/205*100;

A_exp2Algo_001 = find(Expert_label_001(:,2) == Algo_label_001);
A_exp2Algo_001_perc = numel(A_exp2Algo_001)/185*100;

A_exp2Algo_PERC = mean([A_exp2Algo_perc A_exp2Algo_001_perc]);

% AASM vs Algo
A_AASMAlgo = find(CCSHS_label == Algo_label);
A_AASMAlgo_perc = numel(A_AASMAlgo)/205*100;

A_AASMAlgo_001 = find(CCSHS_label_001 == Algo_label_001);
A_AASMAlgo_001_perc = numel(A_AASMAlgo_001)/185*100;

A_AASMAlgo_PERC = mean([A_AASMAlgo_perc A_AASMAlgo_001_perc]);

%%%% Plot
Algo_allExperts = mean([A_exp1Algo_PERC A_exp2Algo_PERC]);
AASM_allExperts = mean([A_exp1AASM_PERC A_exp2AASM_PERC]);

%% std
expexp = [A_expexp_perc A_expexp_001_perc];
expAASM = [A_exp1AASM_perc A_exp1AASM_001_perc A_exp2AASM_perc A_exp2AASM_001_perc];
expFBC = [A_exp1Algo_perc  A_exp1Algo_001_perc A_exp2Algo_perc A_exp2Algo_001_perc];
AASMFBC = [A_AASMAlgo_perc A_AASMAlgo_001_perc];

std_expexp = std(expexp);
std_expAASM = std(expAASM);
std_expFBC = std(expFBC);
std_AASMFBC = std(AASMFBC);

%% ttest (agreemnt sign lower higher?)

[h,p,ci,stats] = ttest2(expexp,expAASM);   % No sign
[h,p,ci,stats] = ttest2(expFBC,AASMFBC);   % No sign


%%

%%%% add error bars first (to hid lower tails)

% Add error bar std (plot it first to hide lower bar and lower end
ALL = [AASM_allExperts; Algo_allExperts; A_AASMAlgo_PERC];
STD = [std_expAASM;std_expFBC;0];

B = bar([AASM_allExperts Algo_allExperts A_AASMAlgo_PERC]);

f = figure; ax=gca; fH = gcf; colormap(jet(4));

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(ALL);
z = B.XEndPoints;

er = errorbar(z,ALL,STD,'k','linestyle','none','LineWidth',1);

hold on 
L = yline(20,'--','LineWidth',1.3);

text(3.5,23,'chance level','FontSize',14)

hold on
B = bar([AASM_allExperts Algo_allExperts A_AASMAlgo_PERC]);

% Parameters

set(gca,'xtick',[])

Y = ylabel('Agreement (%)')
Y.Position = [-0.81959941179007,39.80467795161245,-1];
ylim([0 80])
ax.YGrid = 'on';

%%% Color bars
addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';
[cb] = cbrewer('qual', 'Set2', 12, 'pchip');
index = [1 3 4];
for i = 1:3
    B.FaceColor = 'flat';
    B.CData(i,:) = cb(index(i),:);
end


% legend
%%% Just for legend
hold on; b = bar([NaN],'FaceColor',cb(1,:)); 
hold on; c = bar(NaN,'FaceColor',cb(3,:));
hold on; d = bar(NaN,'FaceColor',cb(4,:));
legend([b c d],'Expert scorers vs. AASM scorers','Expert scorers vs Feature-based clustering','AASM scorers vs Feature-based clustering','Location','southoutside')
legend boxoff


f.Position = [75,212,1248,585];
ax.Position = [0.15,0.25,0.35,0.70];
ax.FontSize = 22;
set(ax,'box','off')

set(gcf,'color','white')
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
fpath = '/Users/nico/Documents/HCTSA/Analysis/ccshs_exp_algo';
%export_fig([fpath filesep 'Agreement'],'-r 300')
