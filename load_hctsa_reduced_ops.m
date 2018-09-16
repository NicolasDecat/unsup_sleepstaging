function [feat_id, feat] = load_hctsa_reduced_ops(filename, hctsafile)
    all_op = load(hctsafile,'Operations');

    fileID = fopen(filename);
    features = textscan(fileID,'%s %s %s');
    fclose(fileID);

    %% Wanted operation names
    feat_name = features{1,2};

    %% Check operation name, get feat_id
    nn=0;
    for n = 1:length(feat_name)
        op_name = char(feat_name(n));
        for i = 1:length(all_op.Operations)
            name = all_op.Operations(i).Name;
            if strcmp(op_name,name)
                nn=nn+1;
                feat_id(nn) = i;
                feat(nn).id = i; % all_op.Operations(i).ID % Actual operation ID
                feat(nn).name = cellstr(name);
            end
        end
    end
end