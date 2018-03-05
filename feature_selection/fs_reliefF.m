function [features, indexes] = fs_reliefF(row, datamat, answer)
    % Assertion strcmp(char(row.fs_algorithm), 'RELIEFF'
    fs_type = char(row.fs_type);
    
    if (size(datamat, 1) ~= length(answer)) 
       disp('Warning. The height of the datamat must match the length of the answer. Datamat could contains multiple channels.');
       datamat = datamat(1:length(answer), :);
    end
    
    [ranks, weights] = relieff(datamat, answer, 20);
    
    selected_feature_indexes = 1:length(ranks);
    if (strcmp(fs_type, 'Top') == 1)
        selected_feature_indexes = 1: row.fs_count;
    elseif (strcmp(fs_type, 'Random') == 1)
        selected_feature_indexes = randperm(length(ranks), row.fs_count);
    end
    
    selected_features = ranks(selected_feature_indexes);

    features = selected_features;
    indexes = selected_feature_indexes;
end