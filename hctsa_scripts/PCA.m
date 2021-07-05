%% Plot tsne

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

Subs = {'005'}; % {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Get the name of all 7749 features 
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Original Labels
    load(sprintf('ccshs_1800%s_annot.mat',sub),'sleepstage')
    labels = sleepstage;
    W = find(labels == 0);
    N1 = find(labels == 1);
    N2 = find(labels == 2);
    N3 = find(labels == 3);
    R = find(labels == 5);
    
    original_labels = [{W},{N1},{N2},{N3},{R}];
    
    
    % Cluster decisions
    load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_3ch_005')
    cluster_decision = statsOut.predictTest;
    CW = find(cluster_decision == 0);
    C1 = find(cluster_decision == 1);
    C2 = find(cluster_decision == 2);
    C3 = find(cluster_decision == 3);
    CR = find(cluster_decision == 5);
    
    Clusters = [{CW} {C1} {C2} {C3} {CR}];

    load('HCTSA_N.mat','TS_DataMat')

    % EEG only
    TS_DataMat = TS_DataMat(1:length(sleepstage),:);

    % t_SNE
    mappedX = tsne(TS_DataMat,'NumDimensions',3);

%     % Plot results    
%      f = figure; ax=gca;
%      
%      for C = 1:5
% 
%          scatter3(mappedX(original_labels{C},1),mappedX(original_labels{C},2),mappedX(original_labels{C},3),50,Colors{C},'filled','LineWidth',1.3)    
%         hold on
% 
%      end

%        % plot their centroid
%      for C = 1:5
%          
%          Centroid = [mean(mappedX(original_labels{C},1)),mean(mappedX(original_labels{C},2)),mean(mappedX(original_labels{C},3))];
%          hold on
%          scatter3(Centroid(1), Centroid(2),Centroid(3),200,'Marker','s','MarkerFaceColor',Colors{C},'MarkerEdgeColor','k','LineWidth',2) 
%          
%      end


     % Plot tSNE
        
     f = figure; ax=gca;
     for C = 1:5

        hold on
        scatter(mappedX(original_labels{C},1),mappedX(original_labels{C},2),50,Colors{C},'filled','LineWidth',1.3)    
    
     end
     
     % plot their centroid
     for C = 1:5
         
         Centroid = [mean(mappedX(Clusters{C},1)) mean(mappedX(Clusters{C},2))];
         hold on
         scatter(Centroid(1), Centroid(2),350,'Marker','s','MarkerFaceColor',Colors{C},'MarkerEdgeColor','k','LineWidth',2) 
         
     end
     
    
    set(gca,'visible','on')
    xlabel('Dimension 1')
    ylabel('Dimension 2')
    ax.FontSize = 18;
    f.Position = [440,78,896,719];
    
       
    %%% Save
    set(f, 'Color', 'w')
    fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
    %export_fig([fpath filesep 'paper_fig_005'],'-r 300')
%     
   
    % Legend
    figure;
    for C = 1:5
    scatter(NaN,NaN,50,Colors{C},'filled','LineWidth',1.3)   
    hold on
    end
    scatter(NaN,NaN,150,'Marker','s','MarkerFaceColor','w','MarkerEdgeColor','k','LineWidth',1) 

    
    legend({'Wake','N1','N2','N3','REM','Centroid'},'FontSize',14);
    legend boxoff

    
end

%% PCA: marker and color

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


Subs = {'005'}; % {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

datasets = [1, 334, 1374; 5, 380, 1442; 439, 150, 1164; 458, 277, 1374; ...
  596, 266, 1396; 604, 502, 1457; 748, 488, 1191; 749,  69, 1009; ...
  752, 147, 1096; 807, 329, 1232; 821, 314, 1316; 870, 282, 1286];


for D = 1:length(Subs)   
    
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
    

    
    load('HCTSA_N.mat','TS_DataMat')

    % EEG only
    TS_DataMat = TS_DataMat(TS,:);

    % t_SNE
    mappedX = tsne(TS_DataMat);

    %%% Plot results

    % Index for each cluster (CW,C1,C2,C3,CR) -> Colors
    % load(sprintf('/Users/nico/Documents/HCTSA/Analysis/Accuracy_100/Matrix_accuracy_per_feat/Table_balanced_%s',sub))
    load(sprintf('/Users/nico/Documents/HCTSA/Analysis/hypnograms/allepochs/statsOut_allepochs_3ch(%s)',sub))
    % cluster_decision = table2array(Table(:,3));
    cluster_decision = statsOut.predictTest;
    cluster_decision = cluster_decision(TS);

    CW = find(cluster_decision == 0);
    C1 = find(cluster_decision == 1);
    C2 = find(cluster_decision == 2);
    C3 = find(cluster_decision == 3);
    CR = find(cluster_decision == 5);
        
    Clusters = [{CW} {C1} {C2} {C3} {CR}];
    
   % Index for each conventional stage (W,N1,N2,N3,R) -> markers
    W = find(labels == 0);
    N1 = find(labels == 1);
    N2 = find(labels == 2);
    N3 = find(labels == 3);
    R = find(labels == 5);
    
    Stages = [{W} {N1} {N2} {N3} {R}];
    
    % All possible cases: rows = CW->CR, col = W->R
      
    for C = 1:5
        for S = 1:5 
        TAB{C,S} = intersect(Clusters{C},Stages{S});
        end
    end
    

    % Markers and colors
    Colors = [{wakeC} {N1C} {N2C} {N3C} {remC}];
    Markers = [{'d'} {'o'} {'^'} {'s'} {'p'}];

    f = figure; ax=gca;
    
     for S = 1:5  
        for C = 1:5    % 5 clusters x 5 stages, potentially 25 combinations of markers / colors
        
            hold on
            scatter(mappedX(TAB{C,S},1),mappedX(TAB{C,S},2),100,Colors{C},Markers{S},'LineWidth',1.3)
        
        end       
     end
     
     % plot their centroid
     for C = 1:5

         Centroid = [mean(mappedX(Clusters{C},1)) mean(mappedX(Clusters{C},2))];
         hold on
         scatter(Centroid(1), Centroid(2),200,'Marker','o','MarkerFaceColor',Colors{C},'MarkerEdgeColor','k','LineWidth',15) 

     end
    set(gca,'visible','on')
    xlabel('Dimension 1')
    ylabel('Dimension 2')
    ax.FontSize = 18;
    f.Position = [440,78,896,719];
    
    %%% Save
    set(f, 'Color', 'w')
    fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
  %  export_fig([fpath filesep 'Mark_Col_005'],'-r 300')

    % Separate file: Fake plots for legend
    figure;
    hold on; bar(NaN,'FaceColor',wakeC)
    hold on; bar(NaN,'FaceColor',N1C)
    hold on; bar(NaN,'FaceColor',N2C)
    hold on; bar(NaN,'FaceColor',N3C)
    hold on; bar(NaN,'FaceColor',remC)
    hold on; plot(NaN,'kd','MarkerSize',10,'LineWidth',1.5)
    hold on; plot(NaN,'ko','MarkerSize',10,'LineWidth',1.5)
    hold on; plot(NaN,'k^','MarkerSize',10,'LineWidth',1.5)
    hold on; plot(NaN,'ks','MarkerSize',10,'LineWidth',1.5)
    hold on; plot(NaN,'kp','MarkerSize',10,'LineWidth',1.5)
    hold on; plot(NaN,'ko','MarkerSize',12,'MarkerFaceColor',[0.7 0.7 0.7],'LineWidth',15)

    legend({'CW','C1','C2','C3','CR','Wake','N1','N2','N3','REM','Centroid'},'FontSize',14);
    legend boxoff
   
    
    fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
    % export_fig([fpath filesep 'Mark_Col_legend'],'-r 300')


%     % Save figure
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
%     saveas(gca,fullfile(fpath,sprintf('PCA(%s)',sub)),'jpg')
    
end



%% Representative PCA for analysis pipeline figure


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

Subs = {'596'}; % {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Get the name of all 7749 features 
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Labels
    load(sprintf('ccshs_1800%s_annot.mat',sub),'sleepstage')
    labels = sleepstage;

    load('HCTSA_N.mat','TS_DataMat')

    % EEG only
    load('/Users/nico/Documents/HCTSA/Analysis/PCA_100/statsOut_3ch_596')
    cluster_decision = statsOut.predictTest;
    TS_DataMat = TS_DataMat(1:length(cluster_decision),:);

    % t_SNE
    mappedX = tsne(TS_DataMat);

    %%% Plot results

    % Index for each cluster (CW,C1,C2,C3,CR) -> Colors
    CW = find(cluster_decision == 0);
    C1 = find(cluster_decision == 1);
    C2 = find(cluster_decision == 2);
    C3 = find(cluster_decision == 3);
    CR = find(cluster_decision == 5);
    
    Clusters = [{CW} {C1} {C2} {C3} {CR}];
    
    %%% Take random values for each cluster (N = 50)

    for C = 1:5
        
        for i = 1:95
             Balanced(i) = randi(length(Clusters{C}));
        end
        
        Clusters{C} = Clusters{C}(Balanced);
    end

    % Markers and colors
    f = figure; ax=gca;
    
     for C = 1:5
         
            hold on

            scatter(mappedX(Clusters{C'},1),mappedX(Clusters{C'},2),80,Colors{C},'Marker','o','MarkerFaceColor',Colors{C})
        
     end
        
end
     
    set(gca,'visible','on')
    ax.FontSize = 18;
    f.Position = [440,78,896,719];
    
    legend({'C1','C2','C3','C4','C5'},'FontSize',16,'Location','southwest');
    legend boxoff
    set(gca,'visible','off')

    %%% Save
    set(f, 'Color', 'w')
    fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
    addpath '/Users/nico/Documents/MATLAB/hctsa-master/export_fig-master'
    export_fig([fpath filesep 'rep_figure_legend'],'-r 300')
    

%% PCA - incorporating time

Subs = {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Get the name of all 7749 features 
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Labels
    load(sprintf('ccshs_1800%s_annot.mat',sub),'sleepstage')
    labels = sleepstage;

    load('HCTSA_N.mat','TS_DataMat')

    %%% Get Time Series from early to late night
    
    % Group TS by stage
    Wake = {find(sleepstage == 0)};
    N1 = {find(sleepstage == 1)};
    N2 = {find(sleepstage == 2)};
    N3 = {find(sleepstage == 3)};
    rem = {find(sleepstage == 5)};
    
    stages = [Wake N1 N2 N3 rem];
    
    % Assign transparency value (from 0.2 to 1) to each data point
    for i = 1:5
        Numshade = length(stages{i});
        Transpar{i} = linspace(1,0.2,Numshade);
    end
    
    % EEG only
    TS_DataMat = TS_DataMat(1:length(sleepstage),:);

    % t_SNE
    mappedX = tsne(TS_DataMat);

    % Plot results
    figure;
    
    ax = axes;
    hold on;

    Colors = [0.894117647058824,0.101960784313725,0.109803921568627;0.215686274509804,0.494117647058824,0.721568627450980;0.301960784313725,0.686274509803922,0.290196078431373;0.596078431372549,0.305882352941177,0.639215686274510;1,0.498039215686275,0];

    for S = 1:5

        Numshade = length(stages{S});   % Number of epochs (of shades)
        whichColor = Colors(S,:);
        Labels = labels(stages{S});

        for P = 1:Numshade
            
            group = Labels(P);
            g = group;
            x = mappedX(stages{S}(P),1);
            y = mappedX(stages{S}(P),2);

            % scatter syntax 
            hold on % important
            uniqueGroups = unique(group);
            h = arrayfun(@(g)scatter(x(group==g), y(group==g), 'filled'), uniqueGroups);
            set(h, 'MarkerFaceColor',whichColor) % set color
            set(h, 'MarkerFaceAlpha', Transpar{S}(P)) % set transparency level

        end

        %%% Create fake graph to adjust legend later
        v(S) = scatter(nan, nan);
        set(v(S), 'MarkerFaceColor',whichColor) 
        set(v(S), 'MarkerEdgeColor','none') 

    end
    
        title(sprintf('Dataset %s',sub))
        legend([v(1) v(2) v(3) v(4) v(5)],{'Wake','N1','N2','N3','REM'});
    
    
%     % Save figure
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100/PCA_time/';
%     exportgraphics(gca,sprintf('/Users/nico/Documents/HCTSA/Analysis/PCA_100/PCA_time/HQ2/PCA(%s).jpg',sub),'Resolution',600)

end


