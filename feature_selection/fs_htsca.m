function [features, indexes] = fs_htsca(row, all_features_count, top_198_features)
    % Assertion strcmp(char(row.fs_algorithm), 'BEN'
    fs_type = char(row.fs_type);
    if (strcmp(fs_type, 'Top') == 1)
        selected_features = top_198_features;
    else
        selected_features = [1:all_features_count];
    end

    selected_feature_indexes = 1:length(selected_features);
    if (strcmp(fs_type, 'Random') == 1)
        selected_feature_indexes = randperm(all_features_count, row.fs_count);
    elseif (strcmp(fs_type, 'Top') == 1)
        selected_feature_indexes = [1: row.fs_count];
    end

    features = selected_features;
    indexes = selected_feature_indexes;
end