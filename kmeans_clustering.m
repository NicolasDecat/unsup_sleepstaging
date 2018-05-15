hctsafile = '180001_HCTSA/HCTSA_N.mat';
all_op = load(hctsafile,'Operations');
OPS_FILE='reduced_ops.txt';
ANSWER_FILE='ccshs_1800001_annot.mat';
%%
annotation = load(ANSWER_FILE);
label = annotation.sleepstage;

fileID = fopen(OPS_FILE);
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
clear i n nn op_name name

%% Use feat_id to select data from full op
datamat = load(hctsafile,'TS_DataMat');
datamat = datamat.TS_DataMat;
datamat = datamat(1:size(datamat,1)/3, :);
datamat = datamat(:,feat_id);

%% Perform clustering
n_clust = 5;
[idx,c, sse] = kmeans(datamat,n_clust,'Distance','sqeuclidean',...
                    'Display','off','Replicates',50,'MaxIter',500);

%% Plot
markersize=12;                
figure;
plot(datamat(idx==1,1),datamat(idx==1,2),'r.','MarkerSize',markersize)
hold on;
plot(datamat(idx==2,1),datamat(idx==2,2),'b.','MarkerSize',markersize)
plot(datamat(idx==3,1),datamat(idx==3,2),'y.','MarkerSize',markersize)
plot(datamat(idx==4,1),datamat(idx==4,2),'g.','MarkerSize',markersize)
plot(datamat(idx==5,1),datamat(idx==5,2),'c.','MarkerSize',markersize)

title("KMeans Clustering");
plot(c(:,1),c(:,2),'kx','MarkerSize',markersize,'LineWidth',2)
%plot(c(:,1),c(:,2),'ko','MarkerSize',12,'LineWidth',2)

legend(['Cluster 1 (SSE: ' num2str(sse(1))],...
       ['Cluster 2 (SSE: ' num2str(sse(2))],...
       ['Cluster 3 (SSE: ' num2str(sse(3))],...
       ['Cluster 4 (SSE: ' num2str(sse(4))],...
       ['Cluster 5 (SSE: ' num2str(sse(5))],...
       'Centroids', 'Location','NW')

label=label+1;
label(label==6)=5;
hold off;

figure;
plot(datamat(label==1,1),datamat(label==1,2),'r.','MarkerSize',markersize)
hold on;
plot(datamat(label==2,1),datamat(label==2,2),'b.','MarkerSize',markersize)
plot(datamat(label==3,1),datamat(label==3,2),'y.','MarkerSize',markersize)
plot(datamat(label==4,1),datamat(label==4,2),'g.','MarkerSize',markersize)
plot(datamat(label==5,1),datamat(label==5,2),'c.','MarkerSize',markersize)
title("Actual Stages");

%plot(c(:,1),c(:,2),'kx','MarkerSize',markersize,'LineWidth',2)
%plot(c(:,1),c(:,2),'ko','MarkerSize',12,'LineWidth',2)

legend(['Cluster 1 (SSE: ' num2str(sse(1))],...
       ['Cluster 2 (SSE: ' num2str(sse(2))],...
       ['Cluster 3 (SSE: ' num2str(sse(3))],...
       ['Cluster 4 (SSE: ' num2str(sse(4))],...
       ['Cluster 5 (SSE: ' num2str(sse(5))],...
       'Centroids', 'Location','NW')

hold off;

