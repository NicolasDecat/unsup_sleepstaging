

%% Accuracy data

% Load data
load('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Per_correct_mean(Dataset 001)')

% Acc per feature, Wake vs N3 classifier  
Data{1,1} = Per_correct_mean(3,:); 
% Acc per feature, Wake vs N1 classifier 
Data{2,1} = Per_correct_mean(1,:); 

%%%% Plot Rainbow graph

% Set colors
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
cl(1, :) = cb(4, :);
cl(2, :) = cb(1, :);
cl(3, :) = cb(2, :);

figure; 
format_fig;

%% Left plot: Wake vs N3

subplot(1,2,2);  

% Wake vs N3: accuracy per feature
h1 = raincloud_plot(Data{1,1}, 'box_on', 1, 'color', cb(7,:), 'alpha', 0.5,...    
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);

yl = ylim;  % used later to get default max on y axis

hold on
i = bar(90.6,yl(2));   % Bar for all-features accuracy (unsup)
hold on 
j = bar(98.2,yl(2));   % Bar for all-features accuracy (sup)
hold on
mean1 = mean(Data{1,1});
a = xline(mean1,'Color',cb(7,:),'LineWidth',2.5);
ylim([0 1])

title('Wake vs N3')
k = line([50 50],[yl(1) yl(2)],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2);  % Line for 50% accuracy
ylim([-0.009 yl(2)]);
legend([h1{1} a i(1) j(1) k(1)], {'Accuracy per feature (unsupervised)','Mean accuracy across features','Accuracy all features (unsup)','Accuracy all features (SVM)','50% accuracy'},'Location','northwest','FontSize',10);


%% Right plot: W vs N1

subplot(1,2,1);  

% Wake vs N1: accuracy per feature
l1 = raincloud_plot(Data{2,1}, 'box_on', 1, 'color', cb(7,:), 'alpha', 0.5,...    
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15, 'box_col_match', 0);

yl = ylim;  % used later to get default max on y axis

hold on
m = bar(70.7,yl(2));   % Bar for all-features accuracy (unsup)
hold on 
n = bar(89.1,yl(2));   % Bar for all-features accuracy (sup)
hold on
mean2 = mean(Data{2,1});
b = xline(mean2,'Color',cb(7,:),'LineWidth',2.5);

title('Wake vs N1')
o = line([50 50],[yl(1) yl(2)],'Color',[1 1 1]*0.7,'LineStyle','--','LineWidth',2);  % Line for 50% accuracy
ylim([-0.015 yl(2)]);
legend([l1{1} b m(1) n(1) o(1)], {'Accuracy per feature (unsupervised)','Mean accuracy across features','Accuracy all features (unsup)','Accuracy all features (SVM)','50% accuracy'},'Location','northwest','FontSize',10);

fpath = '/Users/nico/Documents/HCTSA/Analysis/AUC_100/Figures';
% saveas(gca,fullfile(fpath,sprintf('CEN_all')),'jpeg')


