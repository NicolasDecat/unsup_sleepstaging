
%% Plot average percentage in confusion matrix


% Show or not the figure
set(0,'DefaultFigureVisible','on')

if ~exist('optional_title','var')
      optional_title = "";
end

   scoredTest = statistics.scoredTest;
   predictTest = statistics.predictTest;


%% Confusion matrix of test data
% Reshape scored and predict matrix
g_labelTest = reshape(scoredTest,1,[]);
g_clustTest = reshape(predictTest,1,[]);

% Labelled - make non-zero stage
g_labelTest = g_labelTest+1;

% Clustered - Use final clustering output
g_clustTest = g_clustTest+1;

% Cluster 6 becomes 5
g_labelTest(g_labelTest==6) = 5;
g_clustTest(g_clustTest==6) = 5;

% Visualise confusion matrix
figure;



%% 

answer = g_labelTest;
predict = g_clustTest;

[~, cols] = size(answer);  % length(TS_Mat) * 5

clear max
max_elem = max(answer);    % 5
confmat = zeros(max_elem, max_elem);

for i = 1:cols
    answer_index = answer(:,i); 
    predict_index = predict(:,i);
    confmat(answer_index, predict_index) = confmat(answer_index, predict_index) + 1;
end

totalTarget = sum(confmat, 2);
% perResponse = confmat./repmat(totalTarget, 1, max_elem)*100;
perResponse = MEAN_percent_cf;
t = cell(max_elem, max_elem);
for i=1:max_elem
    for j=1:max_elem
        t(i, j) = cellstr([num2str(round(perResponse(i,j), 2)), '%']);
    end
end

figure;
imagesc(perResponse);

x = repmat(1:max_elem,max_elem,1);
y = x';
text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 12, ...
    'FontWeight', 'bold');
ax = gca;
ax.XTick = 1:max_elem;
ax.YTick = 1:max_elem;
ax.XTickLabels = {'W', 'N1', 'N2', 'N3', 'R'};
ax.YTickLabels = {'W', 'N1', 'N2', 'N3', 'R'};
ylabel('Target Class');
xlabel('Output Class');
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

fpath = '/Users/nico/Documents/HCTSA/Analysis/CF_v2(10iter)/Average CF datasets v2';
saveas(gca,fullfile(fpath,'CF_all_allch'),'jpeg')
% saveas(gca,fullfile(fpath,sprintf('CF_%s_%s(%s)', sub, chan, Nit)),'jpeg')

