
%% Plot the 2D t-SNE of N2 epochs including all subjects

% color
addpath '/Users/nico/Documents/MATLAB/cbrewer/cbrewer/cbrewer';


[RE] = cbrewer('div', 'Spectral', 12, 'pchip'); 
   wakeC = RE(2,:);
   
[BR] = cbrewer('div', 'PuOr', 12, 'pchip');
   N1C = BR(3,:);

[BL] = cbrewer('seq', 'Blues', 12, 'pchip');
   N2C = BL(8,:);
   
[PU] = cbrewer('div', 'PRGn', 12, 'pchip');
   N3C = PU(2,:);

[GR] = cbrewer('seq', 'YlGn', 12, 'pchip');
   remC = GR(7,:);
   
Colors = [{wakeC} {N1C} {N2C} {N3C} {remC}];


Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

datasets = [1, 334, 1374; 5, 380, 1442; 439, 150, 1164; 458, 277, 1374; ...
  596, 266, 1396; 604, 502,21457; 748, 488, 1191; 749,  69, 1009; ...
  752, 147, 1096; 807, 329, 1232; 821, 314, 1316; 870, 282, 1286];



for D = 1:length(Subs)   
%for D = 1

    g = figure; ax=gca;
    
    sub = Subs{D};

    % Get the name of all 7749 features 
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Labels
    load(sprintf('ccshs_1800%s_annot.mat',sub),'sleepstage')
    
    % Remove wakefulness
    LengthTS = numel(sleepstage);
    
    A = find(datasets(:,1) == str2num(sub));
    MIN = datasets(A,2);   % start of sleep (remove period wake)
    MAX = LengthTS;
    TS = MIN:MAX;
    sleepstage = sleepstage(TS);
    labels = sleepstage;
    
    % Get N2 epochs ID
    N2_ID = find(labels == 2);
    
    load('HCTSA_N.mat','TS_DataMat')

    % FV of N2
    TS_DataMat = TS_DataMat(TS,:);
    TS_DataMat = TS_DataMat(N2_ID,:);
    
    mappedX = tsne(TS_DataMat);

    
    % Cluster decisions within N2
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/hypnograms/allepochs/statsOut_allepochs_3ch(%s)',sub))
    cluster_decision = statsOut.predictTest;
    cluster_decision = cluster_decision(TS);
    cluster_decision = cluster_decision(N2_ID);

    CW = find(cluster_decision == 0);
    C1 = find(cluster_decision == 1);
    C2 = find(cluster_decision == 2);
    C3 = find(cluster_decision == 3);
    CR = find(cluster_decision == 5);
        
    Clusters = [{CW} {C1} {C2} {C3} {CR}];

    Colors = [{wakeC} {N1C} {N2C} {N3C} {remC}];

    % t_SNE
    for C = 1:5
        hold on
        scatter(mappedX(Clusters{C},1),mappedX(Clusters{C},2),100,Colors{C},'filled','LineWidth',1.3)
    end
    
    for C = 1:5
        Centroid = [mean(mappedX(Clusters{C},1)) mean(mappedX(Clusters{C},2))];
        hold on
        scatter(Centroid(1), Centroid(2),300,'Marker','s','MarkerFaceColor',Colors{C},'MarkerEdgeColor','k','LineWidth',3) 
    end

    
    set(gca,'visible','on')
%     xlabel('Dimension 1')
%     ylabel('Dimension 2')
     xl = xticklabels;     
     for i = 1:length(xl)
         if bitget(abs(str2num(cell2mat(xl(i)))),1) == 1   % If odd number, remove
             xl{i} = '';
         end
     end

    xticklabels(xl)
    set(gca,'XTickLabel',xticklabels)
    
    yl = yticklabels;
     for i = 1:length(yl)
         if bitget(abs(str2num(cell2mat(yl(i)))),1) == 1   % If odd number, remove
             yl{i} = '';
         end
     end
     
     yticklabels(yl)
     set(gca,'YTickLabel',yticklabels)
     
     xlim([min(mappedX(:,1))-2 max(mappedX(:,1))+2])
     ylim([min(mappedX(:,2))-2 max(mappedX(:,2))+2])


    title(sprintf('%s',sub))

    ax.FontSize = 35;
    g.Position = [389,295,659,502];
    
    %%% Save
    set(g, 'Color', 'w')
    fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
  %  export_fig([fpath filesep 'Mark_Col_005'],'-r 300')

    % Separate file: Fake plots for legend
 
%     hold on; a = scatter(NaN,NaN,'MarkerEdgeColor',wakeC,'LineWidth',15);
%     hold on; b = scatter(NaN,NaN,'MarkerEdgeColor',N1C,'LineWidth',15);
%     hold on; c = scatter(NaN,NaN,'MarkerEdgeColor',N2C,'LineWidth',15);
%     hold on; d = scatter(NaN,NaN,'MarkerEdgeColor',N3C,'LineWidth',15);
%     hold on; e = scatter(NaN,NaN,'MarkerEdgeColor',remC,'LineWidth',15);
%     hold on; f = plot(NaN,'ko','MarkerSize',12,'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',15);
% 
%     title(sprintf('t-SNE on N2 epochs, dataset %s',sub))
%     legend([a b c d e f],{'CW','C1','C2','C3','CR','Centroid'},'FontSize',14);
%     legend boxoff
   
    
    fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100/tSNE_N2';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
    % export_fig([fpath filesep sprintf('%s',sub)],'-r 300')

    
end




%%  SY, SC, SP



