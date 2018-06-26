clc;

basedir = '/Volumes/Spaceship/SleepResearch/sleep_documentation_results/results';
all_dataset = [ ...
  % Male = 1 (First 5) Male = 0 (Rest)
  "1800439", ...
  "1800749", ...
  "1800752", ...
  "1800807", ...
  "1800870", ...
  "1800458", ...
  "1800596", ...
  "1800604", ...
  "1800748", ...
  "1800821", ...
 ];

table_columns = {'Dataset', 'Type', 'NumberOfChannels', 'Accuracy'};
table = array2table(zeros(0,length(table_columns)));
table.Properties.VariableNames = table_columns;

balanced_table_columns = {'Dataset', 'Type', 'NumberOfChannels', 'Accuracy'};
balanced_table = array2table(zeros(0,length(balanced_table_columns)));
balanced_table.Properties.VariableNames = balanced_table_columns;

for i = 1:length(all_dataset)
    dataset_num = str2num(char(strrep(all_dataset(i), '1800', '')));
    exp_result = load(strcat(basedir, filesep, 'sleep_org_', all_dataset(i), '_HCTSA' , filesep, 'SUP_UNSUP_HCTSA_200_DS', ...
        num2str(dataset_num), '.mat'));
    save_stats = exp_result.save_stats;
    balanced_save_stats = save_stats(contains(save_stats.Type, '_BALANCED_LABEL'), :);
    
    [UA, ~, idx] = unique(save_stats(:,[1 6]));
    accum_save_stats = [UA,array2table(accumarray(idx,double(table2array(save_stats(:,4))),[],@mean))];
    row = [ones(size(accum_save_stats, 1), 1) * dataset_num, table2array(accum_save_stats(:, 1)), table2array(accum_save_stats(:, 2)), table2array(accum_save_stats(:, 3))];
    table = [table; array2table(row, 'VariableNames', table_columns)];
    
    balanced_row = [ones(size(balanced_save_stats, 1), 1) * dataset_num, ...
        table2array(balanced_save_stats(:, 1)), ...
        table2array(balanced_save_stats(:, 6)), ...
        table2array(balanced_save_stats(:, 4))];
    
    balanced_table = [balanced_table; array2table(balanced_row, 'VariableNames', balanced_table_columns)];
end

[UA, ~, idx] = unique(table(:,[2 3]));
accum_table = [UA,array2table(accumarray(idx,double(table2array(table(:,4))),[],@mean))];

%% Construct error bar
standard_errors = [];
mean_of_errors = [];

% Filter the type that we want
UA = UA(contains(UA.Type, '_BALANCED_LABEL'), :);

for i = 1:size(UA, 1)
    t = table(table.Type == table2array(UA(i, 1)) & table.NumberOfChannels == table2array(UA(i, 2)), :);
    accuracy = str2double(t.Accuracy);
    disp(accuracy);
    standard_errors = [standard_errors, std(accuracy)./sqrt(size(all_dataset, 2))];
    mean_of_errors = [mean_of_errors; mean(accuracy)];
end


%% Plot 1 - Performance of supervised vs unsupervised (with 6 configurations)
unique_type = unique(accum_table(:, 1));

% Filter the type that we want
unique_type = unique_type(contains(unique_type.Type, '_BALANCED_LABELED'), :);

figure('pos',[250 250 1200 800]);

bar_data = [];
for i = 1:size(unique_type, 1)
   filter = accum_table(accum_table.Type == table2array(unique_type(i, :)), :);

   bar_data = [bar_data; table2array(filter(:, 3))'];
end

bar_data = bar_data.*100;
standard_errors = standard_errors.*100;

set(0,'DefaultTextInterpreter','none');
standard_error_reshape = reshape(standard_errors, 3, 2)';
barplot = barwitherr(standard_error_reshape, bar_data);
hold on;

% barplot = bar(bar_data);
l = cell(1,3);
l{1} = 'EEG';
l{2} = 'EEG + EOG';
l{3} = 'EEG + EOG + EMG';
legend(barplot, l);

title(strcat('Accuracy between supervised and unsupervised algorithms with different combination of channels (n=', num2str(size(all_dataset, 2)), ')'), 'FontSize', 16);
grid on;
xlabel('Algorithm', 'FontSize', 15);
ylabel('Accuracy (%)', 'FontSize', 15);
set(gca, 'xticklabel', [ ...
 "Supervised (Balanced + Labeled)", ...
%  "Supervised (Unbalanced + Labeled A)", ...
%  "Supervised (Unbalanced + Labeled B)", ...
 "Unsupervised (Balanced + Labeled)", ...
%  "Unsupervised (Unbalanced + Labeled A)", ...
%  "Unsupervised (Unbalanced + Labeled B)" 
]);
fix_xticklabels(gca, 0.1, {'FontSize', 14});
set(gca,'TickLabelInterpreter','none');
% xtickangle(25);

% Chance level line
XL = xlim();
plot( XL, [20, 20], '--', 'color', [0.2 0.2 0.2], 'HandleVisibility', 'off');
hold off;

text(0.25, 20, 'Chance level', 'FontSize', 14);

