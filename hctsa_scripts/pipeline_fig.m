
%% Paper figure: analysis pipeline 

%%%%% Analysis pipeline figure (Fig2): desccription of balanced dataset,
%%%%% k-means clustering, cluster mapping and cross-validation


%% Hypnogram + colorbar

% Which dataset?
sub = '596';
cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))

%%%%% Plot the hypnogram

load(sprintf('/Users/nico/Documents/HCTSA/Analysis/PCA_100/statsOut_3ch_%s',sub))
original_labels = statsOut.scoredTest;

% Cut first 2 hours
original_labels = original_labels(240:end);   % because we skipped first 2hrs

% Get idx for each stage. Then give each stage a number, which corresponds to the height of horizontal
% line for the hypnogram
stages = [0 1 2 3 5];
hLine = [180:-20:100];

for S = 1:5
    stage_idx = find(original_labels == stages(S));
    Line_hypno(stage_idx) = hLine(S);
end

NumTS = 1:length(original_labels);

% Plot rectangle
f = figure; ax = gca;
% rectangle('Position', [0 min(hLine)-5 numel(NumTS) max(hLine)-min(hLine)+10], 'FaceColor',[.96 .96 .96],'EdgeColor',[.96 .96 .96]);

% Plot lines
hold on
for L = 1:5
    line([0 numel(NumTS)], [hLine(L) hLine(L)],'LineStyle','-','Color',[.75 .75 .75],'LineWidth',.3);
end

% Plot hypnogram
y = plot(NumTS,Line_hypno,'k','LineWidth',2,'Color',[.3 .3 .3]); 

% y axis
ax.YTick = 100:20:180;
ax.YTickLabels = {'REM','N3','N2','N1','Wake'};
ylim([min(hLine)-25 max(hLine)+5])
xlim([0 numel(NumTS)])

a = get(gca,'YTickLabel');  
set(gca,'YTickLabel',a,'fontsize',14)
f.Position = [47,375,1265,363];

%%%%% Colorbar

% Get the rectangle positions: get length of each bit of stage
[~,Consec] = find(diff(original_labels));
Consec = [Consec numel(original_labels)];   % length of each rectangle
Consec_label = original_labels(Consec);     % Get label of each rectangle
Consec_label_equi = changem(Consec_label,[1 2 3 4 5],[0 1 2 3 5]);

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

for i = 1:length(Consec)

    if i == 1
        rectangle('Position',[0 75 Consec(i) 20],'FaceColor',Colors{Consec_label_equi(i)},'EdgeColor','none');
    else
        rectangle('Position',[Consec(i-1) 75 Consec(i)-Consec(i-1) 20],'FaceColor',Colors{Consec_label_equi(i)},'EdgeColor','none');
    end

end

set(gca,'xtick',[])
oa.YAxisLine = 'off';

addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
set(f, 'Color', 'w')
fpath = '/Users/nico/Documents/HCTSA/Analysis/pipeline_fig';
% export_fig([fpath filesep 'hypno_colorbar'],'-r 300')



%% Shorter, balanced colorbar

g = figure; ax=gca;

Consec = [40:40:200];
Consec_70 = [30:30:160];
Consec_30 = [10:10:50];


for i = 1:length(Consec)

    if i == 1   % 3 lines = 3 series of rectangle (balanced, 70%, 30%
        rectangle('Position',[0 160 Consec(i) 20],'FaceColor',Colors{i},'EdgeColor','none');
        rectangle('Position',[0 130 Consec_70(i) 20],'FaceColor',Colors{i},'EdgeColor','none');
        rectangle('Position',[0 90 Consec_30(i) 20],'FaceColor',Colors{i},'EdgeColor','none');
    else
        rectangle('Position',[Consec(i-1) 160 Consec(i)-Consec(i-1) 20],'FaceColor',Colors{i},'EdgeColor','none');
        rectangle('Position',[Consec_70(i-1) 130 Consec_70(i)-Consec_70(i-1) 20],'FaceColor',Colors{i},'EdgeColor','none');
        rectangle('Position',[Consec_30(i-1) 90 Consec_30(i)-Consec_30(i-1) 20],'FaceColor',Colors{i},'EdgeColor','none');

    end

end

ylim([min(hLine)-25 max(hLine)+5])
xlim([0 numel(NumTS)])

axis off;
g.Position = [47,375,1265,363];

% Legend

hold on; a = bar(NaN,NaN,'FaceColor',wakeC,'EdgeColor','none');
hold on; b = bar(NaN,NaN,'FaceColor',N1C,'EdgeColor','none');
hold on; c = bar(NaN,NaN,'FaceColor',N2C,'EdgeColor','none');
hold on; d = bar(NaN,NaN,'FaceColor',N3C,'EdgeColor','none');
hold on; e = bar(NaN,NaN,'FaceColor',remC,'EdgeColor','none');
legend([a b c d e],'W','N1','N2','N3','R','Location','southoutside','Orientation','horizontal');
legend boxoff
ax.FontSize = 14;

addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
set(g, 'Color', 'w')
fpath = '/Users/nico/Documents/HCTSA/Analysis/pipeline_fig';
export_fig([fpath filesep 'short_colorbar'],'-r 300')


%% K-means clusterning

% color for each cluster
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

h = figure; ax=gca;

% Location of clusters
range_x = [{[.15 .4]} {[.0 .35]} {[.3 .7]} {[.7 .9]} {[.6 .8]}]; % x-axis range
range_y = [{[.85 .5]} {[.5 .1]} {[.5 .3]} {[.4 .9]} {[.8 .95]}]; % y-axis range

Num = 40;  % Plot 80 TS per cluster

% For each cluster, plot random points (TS) within the x and y range
for i = 1:5
    
    % Generate random points 
    rand_x = range_x{i}(1) + (range_x{i}(2)-range_x{i}(1)) .* rand(Num,1);
    rand_y = range_y{i}(1) + (range_y{i}(2)-range_y{i}(1)) .* rand(Num,1);

    scatter(rand_x,rand_y,100,Colors{i},'filled','LineWidth',1.3)
    
    % Get the centroid
    Centroid(i,:) = [mean(rand_x) mean(rand_y)];
  
    hold on;
    
end

% Now, plot some wrong color points

Num = 1;  % Plot 5 TS per cluster
Stage = [1 2 3 4 5];

for Iter = 1:8   % For 5 iterations (5 dots)
    
    for i = 1:5   % Plot one dot somewhere within cluster, random color

        % Generate random points 
        rand_x = range_x{i}(1) + (range_x{i}(2)-range_x{i}(1)) .* rand(Num,1);
        rand_y = range_y{i}(1) + (range_y{i}(2)-range_y{i}(1)) .* rand(Num,1);

        scatter(rand_x,rand_y,100,Colors{Stage(randi(5,1,1))},'filled','LineWidth',1.3) 
        hold on;
        
    end
    
end

% Eventually, plot centroids
for i = 1:5
    scatter(Centroid(i,1), Centroid(i,2),500,'Marker','o','MarkerFaceColor',Colors{i},'MarkerEdgeColor','k','LineWidth',5) 
    hold on; 
end
    
% Plot parameters
xlim([0 1]); ylim([0 1])
axis off;
h.Position = [440,375,563,422];

% Legend, plot fake
a = scatter(NaN, NaN,400,wakeC,'filled','LineWidth',1.3); 
b = scatter(NaN, NaN,400,N1C,'filled','LineWidth',1.3); 
c = scatter(NaN, NaN,400,N2C,'filled','LineWidth',1.3); 
d = scatter(NaN, NaN,400,N3C,'filled','LineWidth',1.3); 
e = scatter(NaN, NaN,400,remC,'filled','LineWidth',1.3); 
f = scatter(NaN, NaN,800,'Marker','o','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','k','LineWidth',1); 
legend([a b c d e f],'W','N1','N2','N3','R','Centroid','Location','eastoutside')
legend boxoff
ax.FontSize = 16;

% Save
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
set(h, 'Color', 'w')
fpath = '/Users/nico/Documents/HCTSA/Analysis/pipeline_fig';
% export_fig([fpath filesep 'legend'],'-r 300')
% save([fpath filesep 'centroids'],'Centroid')




%% Sequential mapping

j = figure;

% Plot centroids
for i = 1:5
    scatter(Centroid(i,1), Centroid(i,2),500,'Marker','o','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','k','LineWidth',5) 
    hold on; 
end
    
% Plot parameters
xlim([0 1]); ylim([0 1])
axis off;
j.Position = [440,375,563,422];

% Save
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
set(j, 'Color', 'w')
fpath = '/Users/nico/Documents/HCTSA/Analysis/pipeline_fig';
% export_fig([fpath filesep 'Centroids_black'],'-r 300')

%% Testing epochs

l = figure; 
clear rand_x rand_y

Num = 10;  % Plot TS per cluster

%%%%%%%%% WITH CENTROID

% For each cluster, plot random points (TS) within the x and y range
for i = 1:5
    
    % Generate random points 
    rand_x(:,i) = range_x{i}(1) + (range_x{i}(2)-range_x{i}(1)) .* rand(Num,1);
    rand_y(:,i) = range_y{i}(1) + (range_y{i}(2)-range_y{i}(1)) .* rand(Num,1);

    scatter(rand_x(:,i),rand_y(:,i),100,Colors{i},'filled','LineWidth',1.3)
  
    hold on;
    
end

% Now, plot some wrong color points

Num = 1;  % Plot 5 TS per cluster
Stage = [1 2 3 4 5];

for Iter = 1:3   % For 5 iterations (5 dots)
    
    for i = 1:5   % Plot one dot somewhere within cluster, random color

        % Generate random points 
        rand_x_W(:,i) = range_x{i}(1) + (range_x{i}(2)-range_x{i}(1)) .* rand(Num,1);
        rand_y_W(:,i) = range_y{i}(1) + (range_y{i}(2)-range_y{i}(1)) .* rand(Num,1);

        randomC(Iter,i) = randi(5,1,1);
        scatter(rand_x_W(:,i),rand_y_W(:,i),100,Colors{Stage(randomC(Iter,i))},'filled','LineWidth',1.3) 
        hold on;
        
    end
    
end

% Plot centroids
load('/Users/nico/Documents/HCTSA/Analysis/pipeline_fig/centroids')
for i = 1:5
    scatter(Centroid(i,1), Centroid(i,2),500,'Marker','o','MarkerFaceColor',Colors{i},'MarkerEdgeColor','k','LineWidth',5) 
    hold on; 
end

% Plot parameters
xlim([0 1]); ylim([0 1])
axis off;
l.Position = [440,375,563,422];

% Save
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
set(l, 'Color', 'w')
fpath = '/Users/nico/Documents/HCTSA/Analysis/pipeline_fig';
export_fig([fpath filesep 'cross_val_Centr'],'-r 300')

%%%%%%%%% WITHOUT CENTROID

m = figure; 

for Iter = 1:3   % For 5 iterations (5 dots)
    
    for i = 1:5   % Plot one dot somewhere within cluster, random color

        scatter(rand_x(:,i),rand_y(:,i),100,Colors{i},'filled','LineWidth',1.3) 
        hold on;
        
    end
    
end

hold on;

for Iter = 1:3   % For 5 iterations (5 dots)
    
    for i = 1:5   % Plot one dot somewhere within cluster, random color

        scatter(rand_x_W(:,i),rand_y_W(:,i),100,Colors{Stage(randomC(Iter,i))},'filled','LineWidth',1.3) 
        hold on;
        
    end
    
end

% Plot parameters
xlim([0 1]); ylim([0 1])
axis off;
m.Position = [440,375,563,422];

% Save
addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
set(m, 'Color', 'w')
fpath = '/Users/nico/Documents/HCTSA/Analysis/pipeline_fig';
export_fig([fpath filesep 'cross_val_wo_Centr'],'-r 300')
