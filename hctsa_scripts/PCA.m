%% Plot PCA - Supervised

% color
[BL] = cbrewer('seq', 'Blues', 12, 'pchip');
N1C = BL(5,:);
N2C = BL(7,:);
N3C = BL(9,:);
[RE] = cbrewer('div', 'Spectral', 12, 'pchip'); 
wakeC = RE(2,:);
[GR] = cbrewer('seq', 'YlGn', 12, 'pchip');
remC = GR(7,:);

Subs = {'439'}; % {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};

    % Get the name of all 7749 features 
    cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub)) 
    
    % Labels
    load(sprintf('ccshs_1800%s_annot.mat',sub),'sleepstage')
    labels = sleepstage;

    load('HCTSA_N.mat','TS_DataMat')

    % EEG only
    TS_DataMat = TS_DataMat(1:length(sleepstage),:);

    % t_SNE
    mappedX = tsne(TS_DataMat);

    % Plot results
    figure;
    p = gscatter(mappedX(:,1), mappedX(:,2),labels,[wakeC;N1C;N2C;N3C;remC]);

    p(1).MarkerSize = 15; p(2).MarkerSize = 15; p(3).MarkerSize = 15; p(4).MarkerSize = 15; p(5).MarkerSize = 15;
    legend({'Wake','N1','N2','N3','REM'},'FontSize',14);
    xticklabels('')
    yticklabels('')
    set(gca,'visible','off')

%     % Save figure
%     fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
%     saveas(gca,fullfile(fpath,sprintf('PCA(%s)',sub)),'jpg')
    
end

%% Plot PCA - Unsupervised (useless)

% First, run hctsa_allfeatures_spectral (1 sub, 1 chan, 100 Nf) % and
% collect for each iteration: 100 epochs ID per stage, their labels (predict_test) and their TS_DataMat
% vector
load('/Users/nico/Documents/HCTSA/Analysis/spectral/statsOut_2')
load('/Users/nico/Documents/HCTSA/Analysis/spectral/testTS_it_2')
load('/Users/nico/Documents/HCTSA/Analysis/spectral/TestMat')

NumIter = 10;
IDX = NumIter*55;

cluster_decisions = statsOut.scoredTest(1:NumIter,:)'; 
EpochID = testTS_it_2(1:NumIter,:)';

% Reshape
cluster_decisions = reshape(cluster_decisions,[IDX,1]);
EpochID = reshape(EpochID,[IDX,1]);

datamat = TestMat';
datamat = cell2mat(datamat);

Subs = {'001'}; % {'001' '005' '439' '458' '596' '748' '749' '752' '604' '807' '821' '870'};

for D = 1:length(Subs)   
    
    sub = Subs{D};
    
    % Labels
    labels = cluster_decisions;

    % EEG only
    TS_DataMat = datamat;

    % t_SNE
    mappedX = tsne(TS_DataMat);

    % Plot results
    figure;
    p = gscatter(mappedX(:,1), mappedX(:,2),labels,[0.894117647058824,0.101960784313725,0.109803921568627;0.215686274509804,0.494117647058824,0.721568627450980;0.301960784313725,0.686274509803922,0.290196078431373;0.596078431372549,0.305882352941177,0.639215686274510;1,0.498039215686275,0]);

    p(1).MarkerSize = 8; p(2).MarkerSize = 8; p(3).MarkerSize = 8; p(4).MarkerSize = 8; p(5).MarkerSize = 8;
    title(sprintf('Dataset %s',sub))
    legend({'Wake','N1','N2','N3','REM'});
    
    % Save figure
    fpath = '/Users/nico/Documents/HCTSA/Analysis/PCA_100';
    % saveas(gca,fullfile(fpath,sprintf('PCA(%s)',sub)),'jpg')
    
end

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


