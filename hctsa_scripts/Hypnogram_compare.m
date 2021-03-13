
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


%%%% Plot the time series of (mis)classified epochs

load('HCTSA_N.mat','TimeSeries')

% Agreed epochs
wake_wake = 248;
N2_N2 = 343;
N3_N3 = 415;

% Misclassified epochs
wake_N2 = 331;
N2_N1 = 843;
N3_N1 = 1026;

figure;
subplot(3,2,1)
TS = TimeSeries.Data(wake_wake,:);
plot(TS)
ylim([-0.3 0.1])
title('wake, wake')

%%% To display hctsa values
% a = imagesc([TS_DataMat(wake_wake,:)',zeros(6006,1),TS_DataMat(wake_N2,:)']);
% cmap(1,1:3) = [.94 .94 .94];   
% colormap(cmap)    

subplot(3,2,2)
TS = TimeSeries.Data(wake_N2,:);
plot(TS)
ylim([-0.3 0.1])
title('wake, N2')

subplot(3,2,3)
TS = TimeSeries.Data(N2_N2,:);
plot(TS)
ylim([-0.3 0.1])
title('N2, N2')

subplot(3,2,4)
TS = TimeSeries.Data(N2_N1,:);
plot(TS)
ylim([-0.3 0.1])
title('N2, N1')

subplot(3,2,5)
TS = TimeSeries.Data(N3_N3,:);
plot(TS)
ylim([-0.3 0.1])
title('N3, N3')

subplot(3,2,6)
TS = TimeSeries.Data(N3_N1,:);
plot(TS)
ylim([-0.3 0.1])
title('N3, N1')



%% Example of Dataset 439, EEG

set(0,'DefaultFigureVisible','on')

% go to 439 directory
cd('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_439')

%%%%%%% Plot the Data matrix (439)
TS_PlotDataMatrix_edited_2('norm')

%%%%%%% Plot the hypnogram
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs(439)')
original_labels = statsOut.scoredTest;
cluster_decision = statsOut.predictTest;

% Give equivalent label for original labels
a = find(original_labels == 0);
b = find(original_labels == 1);
c = find(original_labels == 2);
d = find(original_labels == 3);
e = find(original_labels == 5);

original_labels(a) = 6000;
original_labels(b) = 7000;
original_labels(c) = 8000;
original_labels(d) = 9000;
original_labels(e) = 10000;

% Give equivalent label for original labels
f = find(cluster_decision == 0);
g = find(cluster_decision == 1);
h = find(cluster_decision == 2);
i = find(cluster_decision == 3);
j = find(cluster_decision == 5);

cluster_decision(f) = 6000;
cluster_decision(g) = 7000;
cluster_decision(h) = 8000;
cluster_decision(i) = 9000;
cluster_decision(j) = 10000;

%%%%%%% Plot the hypnogram: Original Labels
x_time = 1:1166;


% Plot
hold on
ax = gca;

y = plot(x_time,original_labels,'r','LineWidth',2.8,'Color',[1 .5 .5]);  % light red

ax.XTick = 1:120:1166;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:12,'uni',0);

ax.YTick = 6000:1000:10000;
ax.YTickLabels = {'Wake','N1','N2','N3','REM'};
ylim([0 10500])

a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',17)

% Cluster decisions
hold on
z = plot(x_time,cluster_decision,'k','LineWidth',1,'Color',[.3 .3 .3]);  % light grey

% Draw the xlines for misclassified epochs

% Agreed epochs
wake_wake = 86 ;
N2_N2 = 476;
N3_N3 = 710;
% Misclassified epochs
wake_N2 = 129 ;
N2_N1 = 531;
N3_N1 = 763;

hold on
line([wake_wake wake_wake], [6000 10500],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
line([wake_N2 wake_N2], [6000 10500],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
line([N2_N2 N2_N2], [6000 10500],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
line([N2_N1 N2_N1], [6000 10500],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
line([N3_N3 N3_N3], [6000 10500],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
line([N3_N1 N3_N1], [6000 10500],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);

hold on
gap = 18;
str = {'1','2','3','4','5','6'};
text(wake_wake-gap,10250,str(1),'FontSize',18,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])
text(wake_N2-gap,10250,str(2),'FontSize',18,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])
text(N2_N2-gap,10250,str(3),'FontSize',18,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])
text(N2_N1-gap,10250,str(4),'FontSize',18,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])
text(N3_N3-gap,10250,str(5),'FontSize',18,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])
text(N3_N1-gap,10250,str(6),'FontSize',18,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])

legend([y z],'Original labels','Cluster decisions','Location','southeastoutside')



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

str = {'1','2','3','4','5','6'};

load('HCTSA_N.mat','TS_DataMat')

figure;  
[ha, pos] = tight_subplot(3,2,[.045 .1],[.1 .05],[.05 .05]);
ax = gca;

for subplot = 1:6

    axes(ha(subplot)); 
    TS = TimeSeries.Data{Epochs(subplot),:};
    plot(TS)
    ylim([-0.7 0.4])
    title(Stages(subplot))
    xticks = 1:192:3840;
    yticks(-0.7:0.2:0.4);
    xticklabels('')
    yticklabels('')
    
    text(80,0.25,str(subplot),'FontSize',20,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745])

end 


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

figure;  
[ha, pos] = tight_subplot(3,2,[.045 .1],[.1 .05],[.05 .05]);
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



%% Other paper figure (data matrix with color bar below)

%%%%%%% Plot the Data matrix (439)
TS_PlotDataMatrix_edited_TS('norm')  % reorder TS to match cluster decisions

% Reorder epochs based on their stages (as labeled by original labels)
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs(439)')
original_labels = statsOut.scoredTest;

wake_OL = find(original_labels == 0);  
N1_OL = find(original_labels == 1);
N2_OL = find(original_labels == 2);
N3_OL = find(original_labels == 3);
rem_OL = find(original_labels == 5);

stage_ordered = [wake_OL N1_OL N2_OL N3_OL rem_OL];

% Reorder epochs based on their stages (as labeled by cluster decisions)
cluster_decision = statsOut.predictTest;

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
end
  
% Perform cumulative addition to get right rect positions
clust_dec_ordered = reshape(clust_dec_ordered',[1 25]);
PosRect = cumsum(clust_dec_ordered')';
PosRect = reshape(PosRect',[5 5]);
clust_dec_ordered = reshape(clust_dec_ordered,[5 5])';
PosRect = PosRect';

% Draw the rectangles: original labels
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');

hold on
a = rectangle('Position',[0 6000 length(wake_OL) 500],'FaceColor',cb(1,:),'EdgeColor','none');
b = rectangle('Position',[a.Position(3) 6000 length(N1_OL) 500],'FaceColor',cb(2,:),'EdgeColor','none');
c = rectangle('Position',[b.Position(1)+b.Position(3) 6000 length(N2_OL) 500],'FaceColor',cb(3,:),'EdgeColor','none');
d = rectangle('Position',[c.Position(1)+c.Position(3) 6000 length(N3_OL) 500],'FaceColor',cb(4,:),'EdgeColor','none');
e = rectangle('Position',[d.Position(1)+d.Position(3) 6000 length(rem_OL) 500],'FaceColor',cb(5,:),'EdgeColor','none');

% Draw the rectangles: cluster decisions

xline(a.Position(3),'k-','Linewidth',2)
xline(b.Position(1)+b.Position(3),'k-','Linewidth',2)
xline(c.Position(1)+c.Position(3),'k-','Linewidth',2)
xline(d.Position(1)+d.Position(3),'k-','Linewidth',2)

hold on

for i = 1:5

    if i == 1
        rectangle('Position',[0 6500 clust_dec_ordered(i,1) 500],'FaceColor',cb(1,:),'EdgeColor','none');
    else
        rectangle('Position',[PosRect(i-1,5) 6500 clust_dec_ordered(i,1) 500],'FaceColor',cb(1,:),'EdgeColor','none');
    end
    
    rectangle('Position',[PosRect(i,1) 6500 clust_dec_ordered(i,2) 500],'FaceColor',cb(2,:),'EdgeColor','none');
    rectangle('Position',[PosRect(i,2) 6500 clust_dec_ordered(i,3) 500],'FaceColor',cb(3,:),'EdgeColor','none');
    rectangle('Position',[PosRect(i,3) 6500 clust_dec_ordered(i,4) 500],'FaceColor',cb(4,:),'EdgeColor','none');
    rectangle('Position',[PosRect(i,4) 6500 clust_dec_ordered(i,5) 500],'FaceColor',cb(5,:),'EdgeColor','none');

end

ylim([0 7000])

% Legend
aline = line(NaN,NaN,'LineWidth',30,'LineStyle','-','Color',cb(1,:));
bline = line(NaN,NaN,'LineWidth',30,'LineStyle','-','Color',cb(2,:));
cline = line(NaN,NaN,'LineWidth',30,'LineStyle','-','Color',cb(3,:));
dline = line(NaN,NaN,'LineWidth',30,'LineStyle','-','Color',cb(4,:));
eline = line(NaN,NaN,'LineWidth',30,'LineStyle','-','Color',cb(5,:));

legend([aline bline cline dline eline],{'Wake','N1','N2','N3','REM'},'Location','southeastoutside','FontSize',20)
labelOrig = ylabel('Original labels');
labelClust = ylabel('Operations');

ax = gca; 
ax.XTick = 1:120:1166;
ax.XTickLabels = strseq('',0:9)';
ax.YTickLabels = strseq('',0:1000:6000);

title('Dataset 439')
xlabel('Time (hours)')
ylabel('Operations')
ax.FontSize = 16;

%% Same; data matrix with top 40 features only


%%% Example of Dataset 001, EEG+EOG+EMG
      
set(0,'DefaultFigureVisible','on')

%%%%%%% Plot the Data matrix
TS_PlotDataMatrix_edited_top40('norm')

%%%%%%% Plot the hypnogram
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs(439)')
original_labels = statsOut.scoredTest;
cluster_decision = statsOut.predictTest;

% Give equivalent label for original labels
a = find(original_labels == 0);
b = find(original_labels == 1);
c = find(original_labels == 2);
d = find(original_labels == 3);
e = find(original_labels == 5);

original_labels(a) = 45;
original_labels(b) = 55;
original_labels(c) = 65;
original_labels(d) = 75;
original_labels(e) = 85;

% Give equivalent label for original labels
f = find(cluster_decision == 0);
g = find(cluster_decision == 1);
h = find(cluster_decision == 2);
i = find(cluster_decision == 3);
j = find(cluster_decision == 5);

cluster_decision(f) = 45;
cluster_decision(g) = 55;
cluster_decision(h) = 65;
cluster_decision(i) = 75;
cluster_decision(j) = 85;

%%%%%%% Plot the hypnogram: Original Labels
x_time = 1:1166;


% Plot
hold on
ax = gca;

y = plot(x_time,original_labels,'r','LineWidth',2.8,'Color',[1 .5 .5]);  % light red

ax.XTick = 1:120:1166;
ax.XTickLabels = arrayfun(@(a)num2str(a),0:12,'uni',0);

ax.YTick = 45:10:85;
ax.YTickLabels = {'Wake','N1','N2','N3','REM'};
ylim([0 90])

a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',14)

% Cluster decisions
hold on
z = plot(x_time,cluster_decision,'k','LineWidth',1,'Color',[.3 .3 .3]);  % light grey

legend([y z],'Original labels','Cluster decisions','Location','southwest')
