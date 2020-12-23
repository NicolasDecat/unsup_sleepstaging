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
    t = cell(max_elem, max_elem);
    for i=1:max_elem
        for j=1:max_elem
            t(i, j) = cellstr([num2str(round(perResponse(i,j), 2)), '%']);
        end
    end
    
    figure;
    imagesc(perResponse);
    title(ttl);

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
end
