function plot_confusion_matrix(experiment, scoredTrain, predictTrain, scoredTest, predictTest, cm_save_dir, optional_title)
    if ~exist('optional_title','var')
          optional_title = "";
    end
 
    %% Confusion matrix of train data
    % Reshape scored and predict matrix
    g_labelTrain = reshape(scoredTrain,1,[]);
    g_clustTrain = reshape(predictTrain,1,[]);

    % Labelled - make non-zero stage
    g_labelTrain = g_labelTrain+1;

    % Clustered - Use final clustering output
    g_clustTrain = g_clustTrain+1;

    % Cluster 6 becomes 5
    g_labelTrain(g_labelTrain==6) = 5;
    g_clustTrain(g_clustTrain==6) = 5;

    % Visualise confusion matrix
    figure;
    plotconfusion_custom(g_labelTrain, g_clustTrain, strcat('Confusion Matrix - Training', optional_title));
    saveas(gcf, strcat(cm_save_dir, filesep, 'CM_TRN_', experiment, '.png'));

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
    plotconfusion_custom(g_labelTest, g_clustTest, strcat('Confusion Matrix - Testing', optional_title));
    saveas(gcf, strcat(cm_save_dir, filesep, 'CM_TST_', experiment, '.png'));
end