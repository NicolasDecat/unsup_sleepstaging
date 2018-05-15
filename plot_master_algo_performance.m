clear all;

WHICH_DATA = 1;

ap = load('MASTER_SLEEP_UNSUP_VS_SUP.mat', 'save_stats');
apt = ap.save_stats;

channels = unique(apt(:,'NumberOfChannels'));
num_of_features = unique(apt(:, 'NumberOfFeatures'));
num_of_features=sort(double(table2array(num_of_features))','a');

types=unique(apt(:, 'Type'));
types=string(table2array(types))';

figure;

for c = 1:size(types', 1) % types
%% Plot output accuracy
accuracy_train = [];
accuracy_test = [];
errors_test = [];

for l=1:size(num_of_features, 2)
    k = num_of_features(l);
    subt = apt(apt.Type == types(c) & double(apt.NumberOfFeatures) == k & double(apt.NumberOfChannels) == 3, :);
    accuracy_train = [accuracy_train, mean(double(table2array(subt(:, 'TrainingAccuracy')))')];
    
    testing_data = double(table2array(subt(:, 'TestingAccuracy'))');
    accuracy_test = [accuracy_test,  mean(testing_data)];
    errors_test = [errors_test, std(testing_data)];
end

%semilogx(x, accuracy_test);
errorbar(double(num_of_features), accuracy_test, errors_test, 'LineWidth', 3, 'CapSize', 15);
set(gca, 'XScale', 'log');
hold on;

end

title(sprintf('Performance between supervised vs unsupervised for Dataset %d (Testing) - EEG only', WHICH_DATA));
ylabel('Accuracy [0-1]');
xlabel('Number of features');
xlim([0 220]);
ylim([0.25 0.90]);

leg = legend(strrep(types, '_', ' '));
set(leg,'Location','northwest');
hold off;

set(0,'DefaultFigureVisible','on') % Uncomment this to enable the figure displaying


