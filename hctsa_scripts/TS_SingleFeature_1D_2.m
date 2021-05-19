function [dataCell] = TS_SingleFeature(whatData,featID,makeViolin,makeNewFigure,opi)
% TS_SingleFeature  Plot distributions for a single feature given a feature ID
%
%---INPUTS:
% whatData: the data to load in (cf. TS_LoadData)
% featID: the ID of the feature to plot
% makeViolin: makes a violin plot instead of overlapping kernel-smoothed distributions
% makeNewFigure: generates a new figure
% whatStat: can provide an already-computed stat for the feature (otherwise will
%           compute a simple linear classification based metric)

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check Inputs:
if nargin < 3
    makeViolin = false;
end
if nargin < 4
    makeNewFigure = false;
end
if nargin < 5
    whatStat = [];
end
if nargin < 6
    beVocal = true;
end

%-------------------------------------------------------------------------------
% Load data:
[TS_DataMat,TimeSeries,Operations,whatDataSource] = TS_LoadData(whatData);
% Get classLabels:
if ismember('Group',TimeSeries.Properties.VariableNames)
    classLabels = categories(TimeSeries.Group);
else
    error('You must assign groups to data to use TS_SingleFeature. Use TS_LabelGroups.');
end
numClasses = length(classLabels);

TS = size(TS_DataMat,1)/7;
TS_DataMat = TS_DataMat(1:TS,:);
TimeSeries = TimeSeries(1:TS,:);

%-------------------------------------------------------------------------------
load('HCTSA_N.mat','Operations')
% op_ind = find(Operations.ID==featID);
op_ind = featID;

if isempty(op_ind)
    error('Operation with ID %u not found in %s',featID,whatDataSource);
end
if beVocal
    fprintf(1,'[%u] %s (%s)\n',featID,Operations.Name{op_ind},Operations.Keywords{op_ind});
end

%-------------------------------------------------------------------------------
% Plot this stuff:
if makeNewFigure
    f = figure('color','w');
end
hold('on')
ax = gca;
% colors = GiveMeColors(numClasses);
% 
[BL] = cbrewer('seq', 'Blues', 12, 'pchip');
N1C = {BL(5,:)};
N2C = {BL(7,:)};
N3C = {BL(9,:)};
[RE] = cbrewer('div', 'Spectral', 12, 'pchip'); 
wakeC = {RE(2,:)};
[GR] = cbrewer('seq', 'YlGn', 12, 'pchip');
remC = {GR(7,:)};

colors = [wakeC;N1C;N2C;N3C;remC]; 

if makeViolin
    dataCell = cell(numClasses,1);
    for i = 1:numClasses
        dataCell{i} = (TS_DataMat(TimeSeries.Group==classLabels{i},op_ind));
    end

    % Re-order groups by mean (excluding any NaNs, descending):
    meanGroup = cellfun(@nanmean,dataCell);
%     [~,ix] = sort(meanGroup,'descend');

   % Re-order groups: wake, N1, N2, N3, REM   
    ix = [1 5 2 3 4];
   
    extraParams = struct();
    extraParams.theColors = colors(ix);
    extraParams.customOffset = -0.5;
    extraParams.offsetRange = 0.7;
    BF_JitteredParallelScatter(dataCell(ix),1,1,0,extraParams);

    % Adjust appearance:
    ax = gca;
    ax.XLim = [0.5+extraParams.customOffset,numClasses+0.5+extraParams.customOffset];
    ax.XTick = extraParams.customOffset+(1:numClasses);
    ax.XTickLabel = {'W','R','N1','N2','N3'};
    % ax.XTickLabelRotation = 30;
    if opi ==1  || opi == 6
        ylabel('Feature value')
    end
    ax.TickLabelInterpreter = 'none';
    if makeNewFigure
        f.Position(3:4) = [402,159];
    end

    % Annotate rectangles for predicted intervals:
    cfnParams = GiveMeDefaultClassificationParams(TimeSeries,[],false);
    cfnParams.numFolds = 0;
    BF_AnnotateRect(TS_DataMat(:,op_ind),TimeSeries.Group,cfnParams,colors,ax,'left');

    % Trim y-limits (with 2% overreach)
    ax.YLim(1) = min(TS_DataMat(:,op_ind))-0.02*range(TS_DataMat(:,op_ind));
    ax.YLim(2) = max(TS_DataMat(:,op_ind))+0.02*range(TS_DataMat(:,op_ind));

    if meanGroup(1) == max(meanGroup)
        set(gca, 'YDir','reverse')
    end

    
else
    linePlots = cell(numClasses,1);
    for i = 1:numClasses
        featVector = TS_DataMat(TimeSeries.Group==classLabels{i},op_ind);
        [~,~,linePlots{i}] = BF_plot_ks(featVector,colors{i},0,2,20,1);
    end
    % Trim x-limits (with 2% overreach)
    ax.XLim(1) = min(TS_DataMat(:,op_ind))-0.02*range(TS_DataMat(:,op_ind));
    ax.XLim(2) = max(TS_DataMat(:,op_ind))+0.02*range(TS_DataMat(:,op_ind));

    % Add a legend:
    legend([linePlots{:}],classLabels,'interpreter','none','Location','best')
    ylabel('Probability density')

    % Annotate rectangles:
    BF_AnnotateRect('diaglinear',TS_DataMat(:,op_ind),TimeSeries.Group,numClasses,colors,ax,'under');

    % Add x-label:
    xlabel('Output')

    % Adjust position
    if makeNewFigure
        f.Position(3:4) = [405,179];
    end
end


end