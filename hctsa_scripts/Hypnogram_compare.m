
%% Plot paper figure hypnogram comparisons

%%% Example of Dataset 001, EEG+EOG+EMG
      
set(0,'DefaultFigureVisible','on')

%%%%%%% Plot the Data matrix
TS_PlotDataMatrix_edited('norm')

%%%%%%% Plot the hypnogram
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs_3ch')
original_labels = statsOut.scoredTest;
cluster_decision = statsOut.predictTest;

% Give equivalent label for original labels
a = find(original_labels == 0);
b = find(original_labels == 1);
c = find(original_labels == 2);
d = find(original_labels == 3);
e = find(original_labels == 5);

original_labels(a) = 18500;
original_labels(b) = 22500;
original_labels(c) = 26500;
original_labels(d) = 30500;
original_labels(e) = 34500;

% Give equivalent label for original labels
f = find(cluster_decision == 0);
g = find(cluster_decision == 1);
h = find(cluster_decision == 2);
i = find(cluster_decision == 3);
j = find(cluster_decision == 5);

cluster_decision(f) = 18500;
cluster_decision(g) = 22500;
cluster_decision(h) = 26500;
cluster_decision(i) = 30500;
cluster_decision(j) = 34500;

%%%%%%% Plot the hypnogram: Original Labels
x_time = 1:1374;


% Plot
hold on
ax = gca;

y = plot(x_time,original_labels,'r','LineWidth',2.8,'Color',[1 .5 .5]);  % light red

ax.XTick = 1:120:1374;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:12,'uni',0);

ax.YTick = 18500:4000:34500;
ax.YTickLabels = {'Wake','N1','N2','N3','REM'};
ylim([0 35000])

a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',14)

% Cluster decisions
hold on
z = plot(x_time,cluster_decision,'k','LineWidth',1,'Color',[.3 .3 .3]);  % light grey

legend([y z],'Original labels','Cluster decisions','Location','southwest')


%% Example of Dataset 439, EEG

set(0,'DefaultFigureVisible','on')

% go to 439 directory
cd('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_439')

%%%%%%% Plot the Data matrix (439)
TS_PlotDataMatrix_edited_2('norm')

%%%%%%% Plot the hypnogram
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs_3ch(439)')
original_labels = statsOut.scoredTest;
cluster_decision = statsOut.predictTest;

% Give equivalent label for original labels
a = find(original_labels == 0);
b = find(original_labels == 1);
c = find(original_labels == 2);
d = find(original_labels == 3);
e = find(original_labels == 5);

original_labels(a) = 6200;
original_labels(b) = 6700;
original_labels(c) = 7200;
original_labels(d) = 7700;
original_labels(e) = 8200;

% Give equivalent label for original labels
f = find(cluster_decision == 0);
g = find(cluster_decision == 1);
h = find(cluster_decision == 2);
i = find(cluster_decision == 3);
j = find(cluster_decision == 5);

cluster_decision(f) = 6200;
cluster_decision(g) = 6700;
cluster_decision(h) = 7200;
cluster_decision(i) = 7700;
cluster_decision(j) = 8200;

%%%%%%% Plot the hypnogram: Original Labels
x_time = 1:1166;


% Plot
hold on
rectangle('Position', [0 5946 1166 8200], 'FaceColor',[.96 .96 .96]);


hold on
LINE = [6200 6700 7200 7700 8200];
for L = 1:5
    line([0 1166], [LINE(L) LINE(L)],'LineStyle','-','Color',[.75 .75 .75],'LineWidth',.3);
end

ST = {'500','1000','1500','2000','2500','3000','3500','4000','4500','5000','5500'};
lines = [500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500];
for L = 1:11
    if L == 1
        text(-30,lines(L),ST(L),'FontSize',18)
    else
    text(-38,lines(L),ST(L),'FontSize',18)
    end
end

hold on
ax = gca;
y = plot(x_time,original_labels,'LineWidth',2,'Color',[.3 .3 .3]);  % light red

ax.XTick = 1:120:1166;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:12,'uni',0);

ax.YTick = 6200:500:8200;
ax.YTickLabels = {'Wake','N1','N2','N3','REM'};
ylim([0 8400])

a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',17)
set(gca, 'TickLength',[0 0])


% Cluster decisions
% hold on
%z = plot(x_time,cluster_decision,'k','LineWidth',1,'Color','r');  % light grey
           

ax.Position = [0.057,0.09,0.88,0.90];

set(gca, 'fontsize', 24);

% Save
% fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
% export_fig([fpath filesep 'hypno_EEG_cl(439)_2'],'-r 300')



%% Dataset 439, EEG+EOG+EMG

TSS = 140;
set(0,'DefaultFigureVisible','on')
sub = '439';

for C = 1:3  % figure for the EEG, then EOG and EMG
    
    CHAN = C;

    cd('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_439')

    %%%%%%% Plot the Data matrix (439)
    [F] = TS_PlotDataMatrix_edited_3(CHAN,'norm');

    %%%%%%% Plot the hypnogram
    load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs_3ch(439)')
    original_labels = statsOut.scoredTest(1+TSS:end);
    cluster_decision = statsOut.predictTest(1+TSS:end);

    % Give equivalent label for original labels
    a = find(original_labels == 0);
    b = find(original_labels == 1);
    c = find(original_labels == 2);
    d = find(original_labels == 3);
    e = find(original_labels == 5);

    original_labels(a) = 8200;
    original_labels(b) = 7700;
    original_labels(c) = 7200;
    original_labels(d) = 6700;
    original_labels(e) = 6200;

    % Give equivalent label for original labels
    f = find(cluster_decision == 0);
    g = find(cluster_decision == 1);
    h = find(cluster_decision == 2);
    i = find(cluster_decision == 3);
    j = find(cluster_decision == 5);

    cluster_decision(f) = 6200;
    cluster_decision(g) = 6700;
    cluster_decision(h) = 7200;
    cluster_decision(i) = 7700;
    cluster_decision(j) = 8200;

    ax = gca;

    ST = {'1000','2000','3000','4000','5000'};
    lines = [1000 2000 3000 4000 5000];
    for L = 1:5
        text(-38,lines(L),ST(L),'FontSize',18)
    end

    ax.XTick = 1:120:1166-TSS;
    ax.XTickLabels = arrayfun(@(a)num2str(a),0:12,'uni',0);

    ax.YTick = 6200:500:8200;

    ylim([0 5946])
    xlim([0 1166-TSS])  %cut first 2 hours of wake
    
    xlabel('')
    xticklabels('')

    a = get(gca,'YTickLabel');  
    set(gca,'YTickLabel',a,'fontsize',17)
    set(gca, 'TickLength',[0 0])

    ax.Position = [0.07,0.09,0.88,0.90];
    % ax.Position = [0.057,0.3,0.88,0.90];  % for colorbar

    set(gca, 'fontsize', 24);

    F.Position = [1,378,1300,419];

    chan = {'EEG','EOG','EMG'};
    h = text(1190-TSS,3250,chan(CHAN),'FontSize',30);
    set(h,'Rotation',90);
    
    % Save
    fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms/SM';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
    % export_fig([fpath filesep sprintf('hypno(%s)_%s',sub,string(chan(C)))],'-r 300')

end

%%%%%%% Plot the hypnogram: Original Labels

    
    f = figure; ax = gca;
    rectangle('Position', [0 5946 1166 8200], 'FaceColor',[.96 .96 .96]);

    x_time = 1+140:1166;

    hold on
    y = plot(x_time,original_labels,'LineWidth',2,'Color',[.3 .3 .3]);  % light red

    hold on
    LINE = [6200 6700 7200 7700 8200];
    for L = 1:5
        line([0 1166], [LINE(L) LINE(L)],'LineStyle','-','Color',[.75 .75 .75],'LineWidth',.3);
    end
    
     ax.YTickLabels = {'REM','N3','N2','N1','Wake'};
     ax.YTick = 6200:500:8200;
     ylim([6000 8400])
     xlim([140 1166])  %cut first 2 hours of wake
     set(gca, 'fontsize', 20);
     ax.XTick = 1+TSS:120:1166;
     ax.XTickLabels = arrayfun(@(a)num2str(a),0:12,'uni',0);
     xlabel('time (hours)')
     set(gca,'box','off') 
     set(gcf,'color','white')
     f.Position = [1,507,1300,290];
     ax.XRuler.Axle.Visible = 'off'; % ax is axis handle

    % Save
    fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms/SM';
    %export_fig([fpath filesep sprintf('HYPNO(%s)',sub)],'-r 300')



%% Plot the time series of (mis)classified epochs

load('HCTSA_N.mat','TimeSeries')

% Agreed epochs
wake_wake = 86 ;
N2_N2 = 476;
N3_N3 = 710;
 
% Misclassified epochs
wake_N2 = 129 ;
N2_N1 = 531;
N3_N1 = 763;

Epochs = [wake_wake wake_N2 N2_N2 N2_N1 N3_N3 N3_N1];
Stages = {'wake, wake','wake, N2','N2, N2','N2, N1','N3, N3','N3, N1'};

Red = {'wake, ','wake, ','N2, ','N2, ','N3, ','N3, '};
Black = {'wake','N2','N2','N1','N3','N1'};

str = {'1','2','3','4','5','6'};

load('HCTSA_N.mat','TS_DataMat')

a = figure;  
[ha, pos] = tight_subplot(3,2,[.090 .1],[.1 .05],[.05 .05]);
ax = gca;

for subplot = 1:6

    axes(ha(subplot)); 
    TS = TimeSeries.Data{Epochs(subplot),:};
    plot(TS)
    ylim([-0.7 0.4])
    title(['\color{red}' Red{subplot}, '\color{black}' Black{subplot}],'FontSize',15)
    xticks([0:1280:3840]);
    yticks(-0.6:0.4:0.2);
    xticklabels({'0','10','20','30'})
%     yticklabels('')
    ylabel('EEG (mV)')
    xlabel('Time (s)')
    
    text(145,-0.55,str(subplot),'FontSize',20,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])

    
end 

%% See

%%% N2 N2 (466) vs N2 wake (464)
%%% N3 N3 (335) vs N3 N2 (331) / N3 REM (761)

figure;  
[ha, pos] = tight_subplot(2,2,[.090 .1],[.1 .05],[.05 .05]);
ax = gca;

corr = 120;
fal = 450;

axes(ha(1)); plot(TimeSeries.Data{corr,:})
axes(ha(2)); plot(TimeSeries.Data{fal,:})
axes(ha(3)); imagesc(TS_DataMat(corr,:)');
axes(ha(4)); imagesc(TS_DataMat(fal,:)');
numColorMapGrads = 6; 
customColorMap = flipud(BF_GetColorMap('redyellowblue',numColorMapGrads,0));
colormap(customColorMap)
%%
% set(a, 'Color', 'w')
% fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
% export_fig([fpath filesep 'hypno_EEG_cl(439)_TS'],'-r 300')


%%% Hctsa values bars

% Agreed epochs
wake_wake = 86 ;
N2_N2 = 476;
N3_N3 = 710;
 
% Misclassified epochs
wake_N2 = 129 ;
N2_N1 = 531;
N3_N1 = 763;


Epochs = [wake_wake wake_N2 N2_N2 N2_N1 N3_N3 N3_N1];

load('HCTSA_N.mat','TS_DataMat')

b = figure;  
[ha, pos] = tight_subplot(3,2,[.090 .1],[.1 .05],[.05 .05]);
ax = gca;

for subplot = 1:6
    
    axes(ha(subplot)); 
    imagesc(TS_DataMat(Epochs(subplot),:)');
    numColorMapGrads = 6; 
    customColorMap = flipud(BF_GetColorMap('redyellowblue',numColorMapGrads,0));
    colormap(customColorMap)
    set(gca,'xtick',[]); set(gca,'ytick',[])
    set(gca,'visible','off')

end

set(b, 'Color', 'w')
% fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
% export_fig([fpath filesep 'hypno_EEG_cl(439)_values'],'-r 300')



%% Other paper figure (data matrix with color bar below)

addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';
%%%%%%% Plot the Data matrix (439)
[f] = TS_PlotDataMatrix_edited_TS('norm');  % reorder TS to match cluster decisions

% Reorder epochs based on their stages (as labeled by original labels)
% load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Table_balanced_439')
% original_labels = table2array(Table(:,2))';

numTS = 151:1163;  %%%% Only EEG (439), without wakefulness before

load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs_3ch(439)')
original_labels = statsOut.scoredTest;
original_labels = original_labels(numTS);

wake_OL = find(original_labels == 0);  
N1_OL = find(original_labels == 1);
N2_OL = find(original_labels == 2);
N3_OL = find(original_labels == 3);
rem_OL = find(original_labels == 5);

stage_ordered = [wake_OL N1_OL N2_OL N3_OL rem_OL];

% Reorder epochs based on their stages (as labeled by cluster decisions)
% cluster_decision = table2array(Table(:,3))';
cluster_decision = statsOut.predictTest;
cluster_decision = cluster_decision(numTS);

stage_ordered = [{wake_OL} {N1_OL} {N2_OL} {N3_OL} {rem_OL}];

% Get the future rectangle positions for cluster decisions
for i = 1:5  
    
    stage = stage_ordered{1,i};
    clust_dec = cluster_decision(stage);
    
    wake_L = length(find(clust_dec == 0));
    N1_L = length(find(clust_dec == 1));
    N2_L = length(find(clust_dec == 2));
    N3_L = length(find(clust_dec == 3));
    rem_L = length(find(clust_dec == 5));

    clust_dec_ordered(i,:) = [wake_L N1_L N2_L N3_L rem_L];
    
    %%% Gets % 
    Wake_per(i) = round(wake_L/numel(stage_ordered{1,i})*100,1);
    N1_per(i) = round(N1_L/numel(stage_ordered{1,i})*100,1);
    N2_per(i) = round(N2_L/numel(stage_ordered{1,i})*100,1);
    N3_per(i) = round(N3_L/numel(stage_ordered{1,i})*100,1);
    rem_per(i) = round(rem_L/numel(stage_ordered{1,i})*100,1);

end
  
% Perform cumulative addition to get right rect positions
clust_dec_ordered = reshape(clust_dec_ordered',[1 25]);
PosRect = cumsum(clust_dec_ordered')';
PosRect = reshape(PosRect',[5 5]);
clust_dec_ordered = reshape(clust_dec_ordered,[5 5])';
PosRect = PosRect';

% Draw the rectangles: original labels
% color
addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';

[RE] = cbrewer('div', 'Spectral', 12, 'pchip'); 
   wakeC = RE(2,:);
   
[BR] = cbrewer('div', 'PuOr', 12, 'pchip');
   N1C = BR(3,:);

[BL] = cbrewer('seq', 'Blues', 12, 'pchip');
   N2C = BL(8,:);
   
[PU] = cbrewer('div', 'PRGn', 12, 'pchip');
   N3C = PU(2,:);

[GR] = cbrewer('seq', 'YlGn', 12, 'pchip');
   remC = GR(7,:);
   
Colors = [{wakeC} {N1C} {N2C} {N3C} {remC}];

hold on

for i = 1:5

    if i == 1
        rectangle('Position',[0 6000 clust_dec_ordered(i,1) 450],'FaceColor',wakeC,'EdgeColor','none');
    else
        rectangle('Position',[PosRect(i-1,5) 6000 clust_dec_ordered(i,1) 450],'FaceColor',wakeC,'EdgeColor','none');
    end
    
    rectangle('Position',[PosRect(i,1) 6000 clust_dec_ordered(i,2) 450],'FaceColor',N1C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,2) 6000 clust_dec_ordered(i,3) 450],'FaceColor',N2C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,3) 6000 clust_dec_ordered(i,4) 450],'FaceColor',N3C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,4) 6000 clust_dec_ordered(i,5) 450],'FaceColor',remC,'EdgeColor','none');

end

ylim([0 6450])

% Legend
aline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',wakeC);
bline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N1C);
cline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N2C);
dline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N3C);
eline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',remC);

legend([aline bline cline dline eline],{'CW','C1','C2','C3','CR'},'Location','southeastoutside','FontSize',22)
legend boxoff

ax = gca; 
% ax.XTickLabels = '';
% ax.YTickLabels = strseq('',1:1000:5946);
% ax.YTick = [1 1000 2000 3000 4000 5000 5946];
ax.YTickLabels = '';

ST = {'500','1000','1500','2000','2500','3000','3500','4000','4500','5000','5500'};
lines = [500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500];
for L = 1:11
    if L == 1
        text(-45,lines(L),ST(L),'FontSize',18)
    else
    text(-59,lines(L),ST(L),'FontSize',18)
    end
end

set(gca,'XTickLabel',[])
set(ax, 'XTick', []);
set(ax,'TickLength',[0 0])

Y = ylabel('Features');
Y.Position = [-76.10583707504806,2990.838068181818,1];
xlabel('Epochs (ordered by sleep stage)')
ax.FontSize = 24;

ax.Position = [0.09,0.05,0.75,0.96];

%%% Change format
f.Position = [0.5 0.5 1100 900];   % [x y width height]

set(gca, 'fontsize', 24);

%%% Save
fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
%export_fig([fpath filesep 'hypnobars_EEG(439)_balanced'],'-r 300')

confmatrix = false; 

if confmatrix == true

    f = figure;
    C_all = [61.0 26.6 3.1 0.4 8.9;10.7 53.3 7.7 0.7 27.6;5.7 15.3 43.1 20.9 15.0;1.6 1.0 16.4 77.4 3.6;6.6 21.3 10.4 1.6 60.2];
    imagesc(C_all)
    colorbar

    for i = 1:5
        for j = 1:5
            t(i,j) = string(C_all(i,j));
        end
    end

    xline(1.5);
    xline(2.5);
    xline(3.5);
    xline(4.5);
    yline(1.5);
    yline(2.5);
    yline(3.5);
    yline(4.5);
    
    
    x = repmat(1:5,5,1);
    y = x';
    text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 20);
    ax = gca;
    ax.XTick = 1:5;
    ax.YTick = 1:5;
    ax.YTickLabels = {'W', 'N1', 'N2', 'N3', 'R'};
    ax.XTickLabels = {'CW', 'C1', 'C2', 'C3', 'CR'};
    set(ax.XAxis,'fontweight','bold');
    set(ax.YAxis,'fontweight','bold');

    Y = ylabel('AASM labels');
    Y.Position = [-0.05875,2.999997615814208,1];
    X = xlabel('Cluster decisions');
    X.Position = [3.000002384185791,-0.017090650231275,1];
    ax.XAxisLocation = 'top';


    c1=[0.8784 0.9176 0.9569]; % W
    c2=[0.5176 0.6471 0.8314]; % light blue
    c3=[0.2667 0.4157 0.7098]; % dark blue
    
    n1=20;
    n2=20;
    cmap=[linspace(c1(1),c2(1),n1);linspace(c1(2),c2(2),n1);linspace(c1(3),c2(3),n1)];
    cmap(:,end+1:end+n2)=[linspace(c2(1),c3(1),n2);linspace(c2(2),c3(2),n2);linspace(c2(3),c3(3),n2)];
    colormap(cmap')
    cB = colorbar;
    cB.Position = [0.885503297963867,0.050515463917526,0.025236593059937,0.8];
    cB.Label.String = 'Percentage overlap';
    cB.Label.FontSize = 20;

    ax.FontSize = 20;
    set(gca,'TickLength',[0 0])
 
    f.Position = [440 298 649 499];
    ax.Position = [0.13 0.05 0.7 0.8];

    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
    set(gcf,'color','white')
    fpath = '/Users/nico/Documents/HCTSA/Analysis/AUC_100/ConfusionMatrix';
  %  export_fig([fpath filesep 'CF_all_3ch_paper_3'],'-r 300')
end

%% Same: but 3 channel derivations

datasets = [1, 334, 1374; 5, 380, 1442; 439, 150, 1164; 458, 277, 1374; ...
  596, 266, 1396; 604, 502, 1457; 748, 488, 1191; 749,  69, 1009; ...
  752, 147, 1096; 807, 329, 1232; 821, 314, 1316; 870, 282, 1286];

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:12

    sub = Subs{D};
    
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))
    
    load(sprintf('ccshs_1800%s_annot.mat',sub), 'sleepstage')
    LengthTS = numel(sleepstage);

    A = find(datasets(:,1) == str2num(sub));
    MIN = datasets(A,2);   % start of sleep (remove period wake)
    MAX = LengthTS;
    
    numTS = [{MIN:MAX} {MIN+LengthTS:MAX*2} {MIN+(LengthTS*2):MAX*3}]; % EEG or EOG or EMG without wakefulness 
    
    
    for C = 3  % figure for the EEG, then EOG and EMG

        CHAN = C;

        addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';

        %%%%%%% Plot the Data matrix (439)
        [f] = TS_PlotDataMatrix_edited_TS(CHAN,numTS,sub,'norm');  % reorder TS to match cluster decisions

        % Reorder epochs based on their stages (as labeled by original labels)
        % load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Table_balanced_439')
        % original_labels = table2array(Table(:,2))';

        load(sprintf('/Users/nico/Documents/HCTSA/Analysis/hypnograms/allepochs/statsOut_allepochs_3ch(%s)',sub))
        original_labels = [statsOut.scoredTest statsOut.scoredTest statsOut.scoredTest];
        original_labels = original_labels(numTS{C});

        wake_OL = find(original_labels == 0);  
        N1_OL = find(original_labels == 1);
        N2_OL = find(original_labels == 2);
        N3_OL = find(original_labels == 3);
        rem_OL = find(original_labels == 5);

        stage_ordered = [wake_OL N1_OL N2_OL N3_OL rem_OL];

        % Reorder epochs based on their stages (as labeled by cluster decisions)
        % cluster_decision = table2array(Table(:,3))';
        cluster_decision = [statsOut.predictTest statsOut.predictTest statsOut.predictTest];
        cluster_decision = cluster_decision(numTS{C});

        stage_ordered = [{wake_OL} {N1_OL} {N2_OL} {N3_OL} {rem_OL}];

        % Get the future rectangle positions for cluster decisions
        for i = 1:5  

            stage = stage_ordered{1,i};
            clust_dec = cluster_decision(stage);

            wake_L = length(find(clust_dec == 0));
            N1_L = length(find(clust_dec == 1));
            N2_L = length(find(clust_dec == 2));
            N3_L = length(find(clust_dec == 3));
            rem_L = length(find(clust_dec == 5));

            clust_dec_ordered(i,:) = [wake_L N1_L N2_L N3_L rem_L];

            %%% Gets % 
            Wake_per(i) = round(wake_L/numel(stage_ordered{1,i})*100,1);
            N1_per(i) = round(N1_L/numel(stage_ordered{1,i})*100,1);
            N2_per(i) = round(N2_L/numel(stage_ordered{1,i})*100,1);
            N3_per(i) = round(N3_L/numel(stage_ordered{1,i})*100,1);
            rem_per(i) = round(rem_L/numel(stage_ordered{1,i})*100,1);

        end

        % Perform cumulative addition to get right rect positions
        clust_dec_ordered = reshape(clust_dec_ordered',[1 25]);
        PosRect = cumsum(clust_dec_ordered')';
        PosRect = reshape(PosRect',[5 5]);
        clust_dec_ordered = reshape(clust_dec_ordered,[5 5])';
        PosRect = PosRect';

        if CHAN == 3

            % Draw the rectangles: original labels
            % color
            addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';

            [RE] = cbrewer('div', 'Spectral', 12, 'pchip'); 
               wakeC = RE(2,:);

            [BR] = cbrewer('div', 'PuOr', 12, 'pchip');
               N1C = BR(3,:);

            [BL] = cbrewer('seq', 'Blues', 12, 'pchip');
               N2C = BL(8,:);

            [PU] = cbrewer('div', 'PRGn', 12, 'pchip');
               N3C = PU(2,:);

            [GR] = cbrewer('seq', 'YlGn', 12, 'pchip');
               remC = GR(7,:);

            Colors = [{wakeC} {N1C} {N2C} {N3C} {remC}];
            
            load('HCTSA_N.mat', 'op_clust')
            Border = numel(op_clust.ord) +150;

            hold on

            for i = 1:5

                if i == 1
                    rectangle('Position',[0 Border clust_dec_ordered(i,1) 900],'FaceColor',wakeC,'EdgeColor','none');
                else
                    rectangle('Position',[PosRect(i-1,5) Border clust_dec_ordered(i,1) 900],'FaceColor',wakeC,'EdgeColor','none');
                end

                rectangle('Position',[PosRect(i,1) Border clust_dec_ordered(i,2) 900],'FaceColor',N1C,'EdgeColor','none');
                rectangle('Position',[PosRect(i,2) Border clust_dec_ordered(i,3) 900],'FaceColor',N2C,'EdgeColor','none');
                rectangle('Position',[PosRect(i,3) Border clust_dec_ordered(i,4) 900],'FaceColor',N3C,'EdgeColor','none');
                rectangle('Position',[PosRect(i,4) Border clust_dec_ordered(i,5) 900],'FaceColor',remC,'EdgeColor','none');

            end

            ylim([0 Border+800])

            % Legend
            aline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',wakeC);
            bline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N1C);
            cline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N2C);
            dline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N3C);
            eline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',remC);

            legend([aline bline cline dline eline],{'CW','C1','C2','C3','CR'},'Location','southoutside','FontSize',22,'Orientation','Horizontal')
            legend boxoff

        end

        ax = gca; 
        % ax.XTickLabels = '';
        % ax.YTickLabels = strseq('',1:1000:5946);
        % ax.YTick = [1 1000 2000 3000 4000 5000 5946];
        ax.YTickLabels = '';

        load('HCTSA_N.mat', 'op_clust')
        
        if numel(op_clust.ord) > 6000
            
            ST = {'1000','2000','3000','4000','5000','6000'};
            lines = [1000 2000 3000 4000 5000 5900];
            for L = 1:6
                text(-45,lines(L),ST(L),'FontSize',18)
            end
            
        else
            
            ST = {'1000','2000','3000','4000','5000'};
            lines = [1000 2000 3000 4000 5000];
            for L = 1:5
                text(-45,lines(L),ST(L),'FontSize',18)
            end
        
        end

        set(gca,'XTickLabel',[])
        set(ax, 'XTick', []);
        set(ax,'TickLength',[0 0])

        Y = ylabel('Features');
        Y.Position = [-62.974355593566656,2990.838068181818,1];

        if CHAN == 3
             xlabel('Epochs (ordered by sleep stage)')
        else
             xlabel('')
        end

        ax.FontSize = 24;
        ax.Position = [0.09,0.15,0.85,0.8];

        chan = {'EEG','EOG','EMG'};
        h = text(numel(numTS{1})+18,3250,chan(CHAN),'FontSize',30);
        set(h,'Rotation',90);

        %%% Change format
        if CHAN == 3
           f.Position = [1,191,1440,680];
        else
           f.Position = [1,191,1440,589];
        end

        set(gca,'box','off')
        ax.XRuler.Axle.Visible = 'off'; % ax is axis handle
        ax.YRuler.Axle.Visible = 'off'; 
        set(gca, 'fontsize', 24);

        chan = {'EEG','EOG','EMG'};

        %%% Save
        addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
        fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms/SM';
        export_fig([fpath filesep sprintf('hypnobars_(%s)_balanced_%s',sub,string(chan(CHAN)))],'-r 300')

    end

end

%% Same; data matrix with top 40 features only
% 
% 
% %%% Example of Dataset 001, EEG+EOG+EMG
%       
% set(0,'DefaultFigureVisible','on')
% 
% %%%%%%% Plot the Data matrix
% TS_PlotDataMatrix_edited_top40('norm')
% 
% %%%%%%% Plot the hypnogram
% load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs(439)')
% original_labels = statsOut.scoredTest;
% cluster_decision = statsOut.predictTest;
% 
% % Give equivalent label for original labels
% a = find(original_labels == 0);
% b = find(original_labels == 1);
% c = find(original_labels == 2);
% d = find(original_labels == 3);
% e = find(original_labels == 5);
% 
% original_labels(a) = 45;
% original_labels(b) = 55;
% original_labels(c) = 65;
% original_labels(d) = 75;
% original_labels(e) = 85;
% 
% % Give equivalent label for original labels
% f = find(cluster_decision == 0);
% g = find(cluster_decision == 1);
% h = find(cluster_decision == 2);
% i = find(cluster_decision == 3);
% j = find(cluster_decision == 5);
% 
% cluster_decision(f) = 45;
% cluster_decision(g) = 55;
% cluster_decision(h) = 65;
% cluster_decision(i) = 75;
% cluster_decision(j) = 85;
% 
% %%%%%%% Plot the hypnogram: Original Labels
% x_time = 1:1166;
% 
% 
% % Plot
% hold on
% ax = gca;
% 
% y = plot(x_time,original_labels,'r','LineWidth',2.8,'Color',[1 .5 .5]);  % light red
% 
% ax.XTick = 1:120:1166;
% ax.XTickLabels = arrayfun(@(a)num2str(a),0:12,'uni',0);
% 
% ax.YTick = 45:10:85;
% ax.YTickLabels = {'Wake','N1','N2','N3','REM'};
% ylim([0 90])
% 
% a = get(gca,'YTickLabel');  
% set(gca,'YTickLabel',a,'fontsize',14)
% 
% % Cluster decisions
% hold on
% z = plot(x_time,cluster_decision,'k','LineWidth',1,'Color',[.3 .3 .3]);  % light grey
% 
% legend([y z],'Original labels','Cluster decisions','Location','southwest')
