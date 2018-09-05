function [s, s_optk, ch, ch_optk, db, db_optk] = calculate_opt_k(datamat, message_prefix)
    clust = zeros(size(datamat,1), 10);
    for i = 1:10
        clust(:, i) = kmeans(datamat, i,'Distance','sqeuclidean',...
                            'Display','off','Replicates',50,'MaxIter',500);
    end

    sil=evalclusters(datamat, clust, 'silhouette');
    c=evalclusters(datamat, clust, 'CalinskiHarabasz');
    dbo=evalclusters(datamat, clust, 'DaviesBouldin');

    total_epochs = size(datamat, 1);

    s = sil.CriterionValues(sil.OptimalK);
    s_optk = sil.OptimalK;
    ch = c.CriterionValues(c.OptimalK);
    ch_optk = c.OptimalK;
    db = dbo.CriterionValues(dbo.OptimalK);
    db_optk = dbo.OptimalK;
    
    fprintf("%s\n------------\n(Number of epochs=%d) \nSilhouette: %.04f (%d) \nCalinskiHarabasz: %.04f (%d) \nDaviesBouldin: %.04f (%d)\n\n", ...
        message_prefix, ...
        total_epochs, ...
        s, ...
        s_optk, ...
        ch, ...
        ch_optk, ...
        db, ...
        db_optk);
end