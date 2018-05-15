clear all;

WHICH_DATA = 0;

ap = load('MASTER_ALGO_PERFORMANCE_2CHANS.mat', 'algo_performance');
apt = ap.algo_performance;

channels = unique(apt(:,'channel'));
exps = unique(apt(:,'exp'));

if WHICH_DATA == 0
    apt_dataset = apt(:, :);
else
    apt_dataset = apt(apt.dataset == WHICH_DATA, :);
end

exps=table2array(exps)';

x = [1,10, 25, 50, 100, 198];
figure;

for c = 1:length(table2array(channels)') % channels
%% Plot output accuracy
accuracy_train = [];
accuracy_test = [];
errors_test = [];

chans=string(cell2mat(table2array(channels)))';

for l=1:length(exps)
    k = exps(l);
    subt = apt_dataset(apt_dataset.channel == chans(c) & apt_dataset.exp ==k, :);
    accuracy_train = [accuracy_train, mean(table2array(subt(:, 'training_percentage'))')];
    
    testing_data = table2array(subt(:, 'testing_percentage'))';
    accuracy_test = [accuracy_test,  mean(testing_data)];
    errors_test = [errors_test, std(testing_data)];
end

%semilogx(x, accuracy_test);
errorbar(x, accuracy_test, errors_test);
set(gca, 'XScale', 'log');
hold on;

end

title(sprintf('Accuracy Report between channels for Dataset %d (Testing)', WHICH_DATA));
ylabel('Accuracy [0-1]');
xlabel('Number of features');
xlim([0 220]);
%ylim([min_scale_y-0.1 max_scale_y+0.1]);

leg = legend(chans);
set(leg,'Location','southeast');
hold off;

set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying


