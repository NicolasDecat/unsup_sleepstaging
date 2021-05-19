function [perResponse] = plotconfusion_custom(answer, predict, ttl)

%    answer = g_labelTest;
%    predict = g_clustTest;
%    ttl = strcat('Confusion Matrix - Testing', optional_title);
    
    [~, cols] = size(answer);
    
    max_elem = max(answer);
    confmat = zeros(max_elem, max_elem);

    for i = 1:cols
        answer_index = answer(:,i); 
        predict_index = predict(:,i);
        confmat(answer_index, predict_index) = confmat(answer_index, predict_index) + 1;
    end
    
    totalTarget = sum(confmat, 2);
    perResponse = confmat./repmat(totalTarget, 1, max_elem)*100;
    
        % CF paper
%     perResponse = [61.0 26.6 3.1 0.4 8.9;...
%         10.7 53.3 7.7 0.7 27.6;...
%         5.7 15.3 43.1 20.9 15.0;...
%         1.6 1.0 16.4 77.4 3.6;...
%         6.6 21.3 10.4 1.6 60.2];
%     
    t = cell(max_elem, max_elem);
    for i=1:max_elem
        for j=1:max_elem
            t(i, j) = cellstr([num2str(round(perResponse(i,j), 1)), '%']);
        end
    end
    
    figure;
    imagesc(perResponse);
    title('');

    x = repmat(1:max_elem,max_elem,1);
    y = x';
    text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 14);
    ax = gca;
    ax.XTick = 1:max_elem;
    ax.YTick = 1:max_elem;
    ax.XTickLabels = {'Wake', 'N1', 'N2', 'N3', 'REM'};
    ax.YTickLabels = {'Wake', 'N1', 'N2', 'N3', 'REM'};
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
    C = colorbar;
    ylabel(C, '% Accuracy','FontSize',15)
    
    ax.FontSize = 14; 

end
