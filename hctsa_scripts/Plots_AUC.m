
%%%%%%%%% Different plots involves in computing the AUC:
%%%%%%%%% 1) to compare OVA, OVO and OVO manual with boxplot and scatter plot (pointless now)
%%%%%%%%% 2) to calculate correlation coefficients for all 3 approaches (pointless now)
%%%%%%%%% 3) to plot testing accuracy for 200 and 7000 features (bar plot) (pointless now)

%%  Barplot of perf across methods

% AUC OVA

% y axis
EEG_OVA = [0.680 0.673 0.689 0.746 0.698 0.729 0.737 0.659 0.655 0.739 0.672 0.553]';
EEG_OVO = [0.698 0.763 0.670 0.667 0.678 0.658 0.693 0.655 0.760 0.776 0.688 0.603]';
EEG_OVOm = [0.688 0.645 0.676 0.583 0.668 0.636 0.685 0.613 0.677 0.774 0.657 0.576]';

EEG2_OVA = [0.786 0.698 0.725 0.774 0.652 0.773 0.771 0.721 0.697 0.753 0.705 0.617]';
EEG2_OVO = [0.762 0.706 0.665 0.716 0.735 0.683 0.704 0.674 0.754 0.696 0.658 0.656]';
EEG2_OVOm = [0.690 0.681 0.639 0.699 0.685 0.652 0.674 0.669 0.675 0.653 0.656 0.640]';

EEG3_OVA = [0.752 0.682 0.718 0.835 0.780 0.781 0.736 0.721 0.697 0.775 0.695 0.591]';
EEG3_OVO = [0.751 0.729 0.627 0.784 0.833 0.784 0.704 0.734 0.752 0.819 0.739 0.577]';
EEG3_OVOm = [0.682 0.694 0.597 0.799 0.845 0.786 0.686 0.664 0.699 0.836 0.738 0.585]';

DataOVA = [EEG_OVA EEG2_OVA EEG3_OVA];  % first columns of group
DataOVO = [EEG_OVO EEG2_OVO EEG3_OVO];   % second columns of group
DataOVOm = [EEG_OVOm EEG2_OVOm EEG3_OVOm];  % tnird columns of group

Data = {DataOVA DataOVO DataOVOm};

figure;
boxplotGroup(Data,'primaryLabels', {'EEG_OVA','EEG_OVO','EEG_OVOm','EEG2_OVA','EEG2_OVO','EEG2_OVOm','EEG3_OVA','EEG3_OVO','EEG3_OVOm'});
title('across methods');
xlabel('methods / channels');
ylabel('AUC');
xtickangle(45)
set(gcf,'color','white')
set(gca,'Fontsize',10);

%% Scatter plot std

x = [1;1;1;2;2;2;3;3;3];

yOVA = [0.053 0.052	0.062]';
yOVO = [0.051 0.036	0.074]';
yOVOm = [0.053 0.020 0.085]';

y = [yOVA; yOVO; yOVOm];

g(1) = {'one-vs-all'}; g(4) = {'one-vs-all'}; g(7) = {'one-vs-all'}; 
g(2) = {'one-vs-one'}; g(5) = {'one-vs-one'}; g(8) = {'one-vs-one'}; 
g(3) = {'one-vs-one manual'}; g(6) = {'one-vs-one manual'}; g(9) = {'one-vs-one manual'}; 
g = g';

figure; gscatter(x,y,g)
xlim([0.5 3.5]);
ylim([0 0.14]);

% Bar graph std

EEG = [0.053 0.051 0.053];  % OVA OVO OVOm
EEG2 = [0.052 0.036 0.020];
EEG3 = [0.062 0.074 0.085];

yOVA = [0.053 0.052	0.062];
yOVO = [0.051 0.036	0.074];
yOVOm = [0.053 0.020 0.085];

y = [EEG; EEG2; EEG3]
figure; bar(y)

%% correlation coefficient

% Testing accuracy:
Sev_TA_EEG = [46.36 38.18 4.55 0 10.91; 8.18 53.64 9.09 0. 29.09; 5.45 34.55 30 10.91 19.09; 1.82 0 22.73 72.73 2.73; 8.18 28.18 10.91 0.91 51.82]; % 7k, EEG
Two_TA_EEG = [62.27 13.64 13.18 0 5.91; 21.82 48.64 10.45 0 19.09; 22.73 16.36 18.64 20.45 21.82; 2.73 0.45 1.36 93.18 2.27; 17.27 40.45 11.82 0.91 29.55]; % 200, EEG
R = corrcoef(Sev_TA_EEG,Two_TA_EEG); % 0.8266

Sev_TA_EEG2 = [58.18 40 0 0 1.82; 4.55 70.91 5.45 0 19.09 ; 4.55 10.91 63.64 11.82 9.09; 1.82 0 3.64 85.45 9.09; 0.91 32.27 17.27 0 44.55]; % 7k, EEG+EOG
Two_TA_EEG2 = [63.18 29.55 0 0 7.27; 16.82 49.09 7.27 0 26.82; 7.27 8.18 66.64 6.82 9.09; 2.27 0 10 87.73 0; 10 20.45 14.55 0 55]; % 200, EEG+EOG
R = corrcoef(Sev_TA_EEG2,Two_TA_EEG2); % 0.955

Sev_TA_EEG3 = [52.73 39.09 6.36 0 1.82; 3.64 63.64 13.64 0.91 18.18; 6.36 18.18 44.55 12.73 18.18; 3.64 2.73 10.91 75.45 7.27; 7.27 20.91 10 4.55 57.27]; % 7k, EEG+EOG+EMG
Two_TA_EEG3 = [54.55 39.09 3.64 0.91 1.82; 14.09 61.82 6.36 0.91 16.82; 9.09 25.91 29.55 9.09 26.36; 2.73 2.27 12.27 79.09 3.64; 7.27 24.55 8.18 0.91 59.09]; % 200, EEG+EOG+EMG
R = corrcoef(Sev_TA_EEG3,Two_TA_EEG3); % 0.989

% one-vs-all
Sev_OVA_EEG = [0.680 0.673 0.689 0.746 0.698 0.729 0.737 0.659 0.655 0.739 0.672 0.553]; % 7k, EEG
Two_OVA_EEG = [0.7159 0.6828 0.6635 0.7349 0.6754 0.7183 0.7411 0.6808 0.6539 0.762 0.6399 0.5387]; % 200, EEG
R = corrcoef(Sev_OVA_EEG,Two_OVA_EEG) % 0.934

% one-vs-all
Sev_OVO_EEG2 = [0.786 0.698 0.725 0.774 0.652 0.773 0.771 0.721 0.697 0.753 0.705 0.617]; % 7k, EEG
Two_OVO_EEG2 = [0.7301 0.6984 0.7293 0.779 0.6911 0.7707 0.796 0.7269 0.6809 0.7587 0.7072 0.6113]; % 200, EEG
R = corrcoef(Sev_OVO_EEG2,Two_OVO_EEG2) % 0.903

% one-vs-all
Sev_OVO_EEG3 = [0.752 0.682 0.718 0.835 0.780 0.781 0.736 0.721 0.697 0.775 0.695 0.591]; % 7k, EEG
Two_OVO_EEG3 = [0.7449 0.7302 0.7313 0.8143 0.7696 0.8005 0.742 0.7399 0.6699 0.775 0.7014 0.5832]; % 200, EEG
R = corrcoef(Sev_OVO_EEG3,Two_OVO_EEG3) % 0.946

% one-vs-one (multiclass.roc)
Sev_OVO_EEG = [0.698 0.763 0.670 0.667 0.678 0.658 0.693 0.655 0.760 0.776 0.688 0.603]; % 7k, EEG
Two_OVO_EEG = [0.657 0.739 0.657 0.677 0.701 0.693 0.686 0.634 0.784 0.829 0.666 0.609]; % 200, EEG
R = corrcoef(Sev_OVO_EEG,Two_OVO_EEG) % 0.898

% one-vs-one (multiclass.roc)
Sev_OVO_EEG2 = [0.762 0.706 0.665 0.716 0.735 0.683 0.704 0.674 0.754 0.696 0.658 0.656]; % 7k, EEG
Two_OVO_EEG2 = [0.766 0.722 0.633 0.692 0.698 0.717 0.677 0.666 0.795 0.698 0.665 0.650]; % 200, EEG
R = corrcoef(Sev_OVO_EEG2,Two_OVO_EEG2) % 0.849

% one-vs-one (multiclass.roc)
Sev_OVO_EEG3 = [0.751 0.729 0.627 0.784 0.833 0.784 0.704 0.734 0.752 0.819 0.739 0.577]; % 7k, EEG
Two_OVO_EEG3 = [0.773 0.775 0.636 0.768 0.850 0.807 0.682 0.770 0.782 0.828 0.720 0.608]; % 200, EEG
R = corrcoef(Sev_OVO_EEG3,Two_OVO_EEG3) % 0.954

% one-vs-one (manual)
Sev_OVOm_EEG = [0.688 0.645 0.676 0.583 0.668 0.636 0.685 0.613 0.677 0.774 0.657 0.576]; % 7k, EEG
Two_OVOm_EEG = [0.668 0.672 0.655 0.627 0.686 0.635 0.678 0.635 0.670 0.789 0.639 0.594]; % 200, EEG
R = corrcoef(Sev_OVOm_EEG,Two_OVOm_EEG) % 0.918

% one-vs-one (manual)
Sev_OVOm_EEG2 = [0.690 0.681 0.639 0.699 0.685 0.652 0.674 0.669 0.675 0.653 0.656 0.640]; % 7k, EEG
Two_OVOm_EEG2 = [0.764 0.701 0.635 0.704 0.680 0.668 0.693 0.675 0.644 0.688 0.664 0.631]; % 200, EEG
R = corrcoef(Sev_OVOm_EEG2,Two_OVOm_EEG2) % 0.719

% one-vs-one (manual)
Sev_OVOm_EEG3 = [0.682 0.694 0.597 0.799 0.845 0.786 0.686 0.664 0.699 0.836 0.738 0.585]; % 7k, EEG
Two_OVOm_EEG3 = [0.776 0.712 0.650 0.797 0.851 0.805 0.693 0.686 0.707 0.832 0.732 0.559]; % 200, EEG
R = corrcoef(Sev_OVOm_EEG3,Two_OVOm_EEG3) % 0.932


%% Bar graph Testing accuracies

Sev_TA_EEG = [0.514 0.467 0.493 0.576 0.504 0.578 0.563 0.481 0.445 0.606 0.461 0.257]; % 7k, EEG
Two_TA_EEG = [0.522 0.476 0.485 0.595 0.502 0.566 0.569 0.495 0.465 0.614 0.454 0.264]; % 200, EEG

Sev_TA_EEG2 = [0.602 0.526 0.573 0.654 0.517 0.637 0.652 0.574 0.481 0.583 0.535 0.382]; % 7k, EEG+EOG
Two_TA_EEG2 = [0.619 0.536 0.560 0.651 0.498 0.633 0.669 0.579 0.491 0.619 0.516 0.384]; % 200, EEG+EOG

Sev_TA_EEG3 = [0.571 0.523 0.553 0.721 0.645 0.680 0.597 0.578 0.482 0.650 0.517 0.313]; % 7k, EEG+EOG+EMG
Two_TA_EEG3 = [0.618 0.543 0.577 0.712 0.643 0.699 0.597 0.585 0.480 0.643 0.517 0.314]; % 200, EEG+EOG+EMG

% Mean 200 features
MEAN_Two_TA_EEG = mean(Two_TA_EEG);
MEAN_Two_TA_EEG2 = mean(Two_TA_EEG2);
MEAN_Two_TA_EEG3 = mean(Two_TA_EEG3);

% Mean 700 features
MEAN_Sev_TA_EEG = mean(Sev_TA_EEG);
MEAN_Sev_TA_EEG2 = mean(Sev_TA_EEG2);
MEAN_Sev_TA_EEG3 = mean(Sev_TA_EEG3);

% bar graph 
X = categorical({'EEG' 'EEG+EOG' 'EEG+EOG+EMG'});
vals = [MEAN_Two_TA_EEG MEAN_Two_TA_EEG2 MEAN_Two_TA_EEG3; MEAN_Sev_TA_EEG MEAN_Sev_TA_EEG2 MEAN_Sev_TA_EEG3];

h = bar(X,vals);
ylim([0.2 0.7])
lgd = legend(h,'200 features','7000 features');
lgd.FontSize = 12;
title('Testing accuracy')

