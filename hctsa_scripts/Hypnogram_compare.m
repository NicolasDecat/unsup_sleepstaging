
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
        text(-27,lines(L),ST(L),'FontSize',14)
    else
    text(-35,lines(L),ST(L),'FontSize',14)
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
z = plot(x_time,cluster_decision,'k','LineWidth',1,'Color','r');  % light grey

%%% Same figure but 4th-7th hour (magnify)

% a = figure;
% copyobj(ax,a)
% set(a, 'Color', 'w')
% customColorMap = flipud(BF_GetColorMap('redyellowblue',6,0));
% colormap(customColorMap);
% 
% factor = (2*60);  % multiply by 2 to go 30-->60s, and *60 for hour
% 
% xlim([4*factor 7*factor ])  % 4th to 5th hour



% Draw the xlines for misclassified epochs

% Agreed epochs
wake_wake = 86;
N2_N2 = 476;
N3_N3 = 710;
% Misclassified epochs
wake_N2 = 129;
N2_N1 = 531;
N3_N1 = 763;

% hold on
% line([wake_wake wake_wake], [6000 7930],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
% line([wake_N2 wake_N2], [6000 7930],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
% line([N2_N2 N2_N2], [6000 7930],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
% line([N2_N1 N2_N1], [6000 7930],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
% line([N3_N3 N3_N3], [6000 7930],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
% line([N3_N1 N3_N1], [6000 7930],'LineStyle',':','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'LineWidth',3);
% 
% hold on
% gap = 4.5;
% str = {'1','2','3','4','5','6'};
% text(wake_wake-gap,8000,str(1),'FontSize',15,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'BackgroundColor','k')
% text(wake_N2-gap,8000,str(2),'FontSize',15,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'BackgroundColor','k')
% text(N2_N2-gap,8000,str(3),'FontSize',15,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'BackgroundColor','k')
% text(N2_N1-gap,8000,str(4),'FontSize',15,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'BackgroundColor','k')
% text(N3_N3-gap,8000,str(5),'FontSize',15,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'BackgroundColor','k')
% text(N3_N1-gap,8000,str(6),'FontSize',15,'FontWeight','bold','Color',[0.552941176470588,0.827450980392157,0.780392156862745],'BackgroundColor','k')

% legend([y],'\newlineOriginal\newlinelabels','Location','southeastoutside')
% legend boxoff               

ax.Position = [0.048,0.08,0.90,0.90];

% Save
% fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
% export_fig([fpath filesep 'hypno_EEG_cl(439)_2'],'-r 300')



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

%%%%%%% Plot the Data matrix (439)
[f] = TS_PlotDataMatrix_edited_TS('norm');  % reorder TS to match cluster decisions

% Reorder epochs based on their stages (as labeled by original labels)
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs_3ch(439)')
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
[BL] = cbrewer('seq', 'Blues', 12, 'pchip');
N1C = BL(5,:);
N2C = BL(7,:);
N3C = BL(9,:);
[RE] = cbrewer('div', 'Spectral', 12, 'pchip'); 
wakeC = RE(2,:);
[GR] = cbrewer('seq', 'YlGn', 12, 'pchip');
remC = GR(7,:);



hold on
a = rectangle('Position',[0 6000 length(wake_OL) 350],'FaceColor',wakeC,'EdgeColor','none');
b = rectangle('Position',[a.Position(3) 6000 length(N1_OL) 350],'FaceColor',N1C,'EdgeColor','none');
c = rectangle('Position',[b.Position(1)+b.Position(3) 6000 length(N2_OL) 350],'FaceColor',N2C,'EdgeColor','none');
d = rectangle('Position',[c.Position(1)+c.Position(3) 6000 length(N3_OL) 350],'FaceColor',N3C,'EdgeColor','none');
e = rectangle('Position',[d.Position(1)+d.Position(3) 6000 length(rem_OL) 350],'FaceColor',remC,'EdgeColor','none');

% Draw the rectangles: cluster decisions

% xline(a.Position(3),'k-','Linewidth',2)
% xline(b.Position(1)+b.Position(3),'k-','Linewidth',2)
% xline(c.Position(1)+c.Position(3),'k-','Linewidth',2)
% xline(d.Position(1)+d.Position(3),'k-','Linewidth',2)


hold on

for i = 1:5

    if i == 1
        rectangle('Position',[0 6350 clust_dec_ordered(i,1) 350],'FaceColor',wakeC,'EdgeColor','none');
    else
        rectangle('Position',[PosRect(i-1,5) 6350 clust_dec_ordered(i,1) 350],'FaceColor',wakeC,'EdgeColor','none');
    end
    
    rectangle('Position',[PosRect(i,1) 6350 clust_dec_ordered(i,2) 350],'FaceColor',N1C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,2) 6350 clust_dec_ordered(i,3) 350],'FaceColor',N2C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,3) 6350 clust_dec_ordered(i,4) 350],'FaceColor',N3C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,4) 6350 clust_dec_ordered(i,5) 350],'FaceColor',remC,'EdgeColor','none');

end

ylim([0 6700])

% Legend
aline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',wakeC);
bline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N1C);
cline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N2C);
dline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N3C);
eline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',remC);

legend([aline bline cline dline eline],{'Wake','N1','N2','N3','REM'},'Location','southeastoutside','FontSize',10)

text(-128,6185,'Original labels','FontSize',11,'FontWeight','bold','Color',[.3 .3 .3]);
text(-154,6500,'Cluster decisions','FontSize',11,'FontWeight','bold','Color',[.3 .3 .3]);


ax = gca; 
% ax.XTickLabels = '';
% ax.YTickLabels = strseq('',1:1000:5946);
ax.YTick = [1 1000 2000 3000 4000 5000 5946];
ax.YTickLabels = arrayfun(@(a)num2str(a),[1 1000 2000 3000 4000 5000 5946],'uni',0);

set(gca,'XTickLabel',[])
set(ax, 'XTick', []);
set(ax,'TickLength',[0 0])

ylabel('Operations')
xlabel('Epochs')
ax.FontSize = 12;

ax.Position = [0.11,0.03,0.8,0.96];
pos = get(gcf,'Position');

%%% Change format
f.Position = [0.5 0.5 900 900];   % [x y width height]

%%% Save
% fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
% export_fig([fpath filesep 'hypnobars_EEG(439)'],'-r 300')
confmatrix = true; 

if confmatrix == true
    
    figure;
    C = [Wake_per;N1_per;N2_per;N3_per;rem_per]';
    C_all = [53.2 27.1 6.6 0.5 12.7; 9.3 48.6 7.7 1.0 33.5; 5.0 14.1 43.0 18.6 19.3; 1.3 1.4 14.1 78.9 4.4; 5.6 28.5 10.5 1.4 54.0];
    imagesc(C_all)
    colorbar

    for i = 1:5
        for j = 1:5
            t(i,j) = {strcat(num2str(C_all(i,j)),' %')};
        end
    end


    x = repmat(1:5,5,1);
    y = x';
    text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 12, ...
        'FontWeight', 'bold');
    ax = gca;
    ax.XTick = 1:5;
    ax.YTick = 1:5;
    ax.XTickLabels = {'W', 'N1', 'N2', 'N3', 'R'};
    ax.YTickLabels = {'W', 'N1', 'N2', 'N3', 'R'};
    ylabel('Original labels');
    xlabel('Cluster decisions');
    ax.XAxisLocation = 'top';

    %Define colormap
    c1=[0 0.65 0]; %G
    c2=[1 1 0]; %Y
    c3=[1 0 0]; %R
    n1=20;
    n2=20;
    cmap=[linspace(c1(1),c2(1),n1);linspace(c1(2),c2(2),n1);linspace(c1(3),c2(3),n1)];
    cmap(:,end+1:end+n2)=[linspace(c2(1),c3(1),n2);linspace(c2(2),c3(2),n2);linspace(c2(3),c3(3),n2)];
    colormap(cmap')
    colorbar

end

%% Same: but 3 channel derivations


%%%%%%% Plot the Data matrix (439)
[f] = TS_PlotDataMatrix_edited_TS('norm');  % reorder TS to match cluster decisions

% Reorder epochs based on their stages (as labeled by original labels)
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs_3ch(439)')
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
[BL] = cbrewer('seq', 'Blues', 12, 'pchip');
N1C = BL(5,:);
N2C = BL(7,:);
N3C = BL(9,:);
[RE] = cbrewer('div', 'Spectral', 12, 'pchip'); 
wakeC = RE(2,:);
[GR] = cbrewer('seq', 'YlGn', 12, 'pchip');
remC = GR(7,:);


hold on
a = rectangle('Position',[0 18000 length(wake_OL) 800],'FaceColor',wakeC,'EdgeColor','none');
b = rectangle('Position',[a.Position(3) 18000 length(N1_OL) 800],'FaceColor',N1C,'EdgeColor','none');
c = rectangle('Position',[b.Position(1)+b.Position(3) 18000 length(N2_OL) 800],'FaceColor',N2C,'EdgeColor','none');
d = rectangle('Position',[c.Position(1)+c.Position(3) 18000 length(N3_OL) 800],'FaceColor',N3C,'EdgeColor','none');
e = rectangle('Position',[d.Position(1)+d.Position(3) 18000 length(rem_OL) 800],'FaceColor',remC,'EdgeColor','none');

% Draw the rectangles: cluster decisions

% xline(a.Position(3),'k-','Linewidth',2)
% xline(b.Position(1)+b.Position(3),'k-','Linewidth',2)
% xline(c.Position(1)+c.Position(3),'k-','Linewidth',2)
% xline(d.Position(1)+d.Position(3),'k-','Linewidth',2)


hold on

for i = 1:5

    if i == 1
        rectangle('Position',[0 18800 clust_dec_ordered(i,1) 800],'FaceColor',wakeC,'EdgeColor','none');
    else
        rectangle('Position',[PosRect(i-1,5) 18800 clust_dec_ordered(i,1) 800],'FaceColor',wakeC,'EdgeColor','none');
    end
    
    rectangle('Position',[PosRect(i,1) 18800 clust_dec_ordered(i,2) 800],'FaceColor',N1C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,2) 18800 clust_dec_ordered(i,3) 800],'FaceColor',N2C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,3) 18800 clust_dec_ordered(i,4) 800],'FaceColor',N3C,'EdgeColor','none');
    rectangle('Position',[PosRect(i,4) 18800 clust_dec_ordered(i,5) 800],'FaceColor',remC,'EdgeColor','none');

end

ylim([0 19600])

% Legend
aline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',wakeC);
bline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N1C);
cline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N2C);
dline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',N3C);
eline = line(NaN,NaN,'LineWidth',25,'LineStyle','-','Color',remC);

legend([aline bline cline dline eline],{'Wake','N1','N2','N3','REM'},'Location','southeastoutside','FontSize',10)

% text(-128,18300,'Original labels','FontSize',11,'FontWeight','bold','Color',[.3 .3 .3]);
% text(-154,18800,'Cluster decisions','FontSize',11,'FontWeight','bold','Color',[.3 .3 .3]);
text(-190,18400,'Original labels','FontSize',11,'FontWeight','bold','Color',[.3 .3 .3]);
text(-220,19000,'Cluster decisions','FontSize',11,'FontWeight','bold','Color',[.3 .3 .3]);


ax = gca; 
% ax.XTickLabels = '';
% ax.YTickLabels = strseq('',1:1000:5946);
ax.YTick = [1000 3000 5000 5946+1000 5946+3000 5946+5000 5946*2+1000 5946*2+3000 5946*2+5000];
ax.YTickLabels = arrayfun(@(a)num2str(a),[1000 3000 5000 1000 3000 5000 1000 3000 5000],'uni',0);

set(gca,'XTickLabel',[])
set(ax, 'XTick', []);
set(ax,'TickLength',[0 0])

ylabel('Operations')
xlabel('Epochs')
ax.FontSize = 12;

% ax.Position = [0.11,0.09,0.8,0.88];
ax.Position = [0.15,0.03,0.72,0.96];

pos = get(gcf,'Position');

%%% Change format
f.Position = [0.5 0.5 900 900];   % [x y width height]

%%% Save
% fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
% export_fig([fpath filesep 'hypnobars_EEG_3ch(439)'],'-r 300')


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
