
%% Paper figure: analysis pipeline 

%%%%% Analysis pipeline figure (Fig2): desccription of balanced dataset,
%%%%% k-means clustering, cluster mapping and cross-validation


%% Hypnogram + colorbar

%%%%% Hypnogram 

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



