
%%%%%%% Had to change l.114: Change TimeSeries.Data into cell array 

function [F] = TS_PlotDataMatrix_edited_3(CHAN,varargin)

TSS = 140;

% TS_PlotDataMatrix   Plot the data matrix.
%
%---EXAMPLE USAGE:
% TS_PlotDataMatrix; % plots clustered data if it exists
% TS_PlotDataMatrix('whatData','norm'); % plots normalized data
%
%---INPUTS:
% whatData: specify 'norm' for normalized data in HCTSA_N.mat, 'cl' for clustered
%         data in HCTSA_cl.mat (default), or specify a filename to load data
%         from that file.
% addTimeSeries: set to 1 to annotate time series segments to the left of the data matrix
% timeSeriesLength: length of time-series annotations (number of samples)
% colorGroups: set to 1 to color time-series groups with different colormaps
% customColorMap: use a custom color map (name to match an option in BF_GetColorMap)
% colorNaNs: whether to plot rectangles over special-values in the matrix (default: 1)
% customOrder: reorder rows and columns according to provided permutation vectors
%
%---OUTPUT:
% Produces a heat map of the data matrix with time series as rows and
%   operations as columns.

% ------------------------------------------------------------------------------

%% Check inputs and set defaults:
% ------------------------------------------------------------------------------
inputP = inputParser;

% whatDataFile
default_whatData = 'norm';
check_whatData = @(x) ischar(x) || isstruct(x);
addOptional(inputP,'whatData',default_whatData,check_whatData);

% addTimeSeries, annotates time series segments to the side of the plot
default_addTimeSeries = true;
check_addTimeSeries = @(x) (isnumeric(x) || islogical(x)) && (x==0 || x==1);
addOptional(inputP,'addTimeSeries',default_addTimeSeries,check_addTimeSeries);

% timeSeriesLength, length of time-series annotations to the left of the main plot
default_timeSeriesLength = 200;
addOptional(inputP,'timeSeriesLength',default_timeSeriesLength,@isnumeric);

% colorGroups, color groups of time series differently:
default_colorGroups = false;
check_colorGroups = @(x) (x==0 || x==1);
addOptional(inputP,'colorGroups',default_colorGroups,check_colorGroups);

% groupReorder, reorder within groups of time series:
default_groupReorder = false;
check_groupReorder = @(x) (x==0 || x==1);
addOptional(inputP,'groupReorder',default_groupReorder,check_groupReorder);

% custom color map, customColorMap
default_customColorMap = 'redyellowblue';
addOptional(inputP,'customColorMap',default_customColorMap,@ischar);

% colorNaNs
default_colorNaNs = 1;
check_colorNaNs = @(x) (x==0 || x==1);
addOptional(inputP,'colorNaNs',default_colorNaNs,check_colorNaNs);

% customOrder (reorder before plotting)
default_customOrder = {[],[]};
check_customOrder = @(x)iscell(x) && length(x)==2;
addOptional(inputP,'customOrder',default_customOrder,check_customOrder);

%% Parse inputs:
parse(inputP,varargin{:});

% Make variables from results of input parser:
addTimeSeries = inputP.Results.addTimeSeries;
whatData = inputP.Results.whatData;
timeSeriesLength = inputP.Results.timeSeriesLength;
colorGroups = inputP.Results.colorGroups;
groupReorder = inputP.Results.groupReorder;
customColorMap = inputP.Results.customColorMap;
colorNaNs = inputP.Results.colorNaNs;
customOrder = inputP.Results.customOrder;
clear('inputP');

% --------------------------------------------------------------------------
%% Read in the data
% --------------------------------------------------------------------------
% You always want to retrieve and plot the clustered data if it exists
getClustered = false;
[TS_DataMat,TimeSeries,Operations] = TS_LoadData(whatData,getClustered);
[numTS,numOps] = size(TS_DataMat); % size of the data matrix


%%%%%%%% Nico modification: Change TimeSeries.Data into cell array (ONLY
%%%%%%%% FOR HCTSA_N FILES -- HCTSA.mat files already in cell array)

CellData = [];
for D = 1:size(TimeSeries.Data,1)
    CellData{D,1} = TimeSeries.Data(D,:);
end

TimeSeries.Data = CellData;

%%%% EEG, EOG and EMG


if CHAN == 1
    TimeSeries = TimeSeries(1+140:1166,:);
    TS_DataMat = TS_DataMat(1+140:1166,:);
elseif CHAN == 2
    TimeSeries = TimeSeries(1167+140:2332,:);
    TS_DataMat = TS_DataMat(1167+140:2332,:);
elseif CHAN == 3
    TimeSeries = TimeSeries(2333+140:3498,:);
    TS_DataMat = TS_DataMat(2333+140:3498,:);
end
% ------------------------------------------------------------------------------
%% Reorder according to customOrder
% ------------------------------------------------------------------------------
if ~isempty(customOrder{1}) % reorder rows
	fprintf(1,'Reordering time series according to custom order specified.\n');
	TS_DataMat = TS_DataMat(customOrder{1},:);
    TimeSeries = TimeSeries(customOrder{1},:);
end
if ~isempty(customOrder{2}) % reorder columns
	fprintf(1,'Reordering operations according to custom order specified.\n');
	TS_DataMat = TS_DataMat(:,customOrder{2});
    Operations = Operations(customOrder{2},:);
end

%-------------------------------------------------------------------------------
% Check group information
%-------------------------------------------------------------------------------
if ismember('Group',TimeSeries.Properties.VariableNames)
	timeSeriesGroups = TimeSeries.Group;
	numClasses = length(categories(timeSeriesGroups));
else
	timeSeriesGroups = [];
end
if colorGroups
	if ~isempty(timeSeriesGroups)
	    fprintf(1,'Coloring time series by group assignment...\n');
	else
	    warning('No group information found')
	    colorGroups = false;
	end
end

%-------------------------------------------------------------------------------
% Reorder according to groups
%-------------------------------------------------------------------------------
groupReorder=false;

if groupReorder
	if isempty(timeSeriesGroups)
		warning('Cannot reorder by time series group; no group information found')
	else
	    [~,ixData] = sort(timeSeriesGroups);
	    dataMatReOrd = TS_DataMat(ixData,:);
	    ixAgain = ixData;
	    for i = 1:numClasses
	        isGroup = grp2idx(TimeSeries.Group(ixData))==i;
	        ordering = BF_ClusterReorder(dataMatReOrd(isGroup,:),'euclidean','average');
	        istmp = ixData(isGroup);
	        ixAgain(isGroup) = istmp(ordering);
	    end
	    ixData = ixAgain; % set ordering to ordering within groups
	    TimeSeries = TimeSeries(ixData,:);
	    TS_DataMat = TS_DataMat(ixData,:);
	    timeSeriesGroups = TimeSeries.Group;
	end
end

% --------------------------------------------------------------------------
%% Prepare data matrix for plotting
% --------------------------------------------------------------------------
if colorGroups
    gi = BF_ToGroup(timeSeriesGroups);

    numGroups = length(gi);

    % Add a group for unlabelled data items if they exist
    if sum(cellfun(@length,gi)) < numTS
        % Add an unlabelled class
        gi0 = gi;
        gi = cell(numGroups+1,1);
        for i = 1:numGroups
            gi{i} = gi0{i};
        end
        clear gi0;
        gi{end} = setxor(1:numTS,vertcat(gi{:}));
        numGroups = numGroups + 1;
    end

    fprintf(1,'Coloring data according to %u groups\n',numGroups);

    % Change range of TS_DataMat to make use of new colormap appropriately
    almost1 = 1 - 1e-7;
    squashMe = @(x) almost1*(x - min(x(:)))/(max(x(:))-min(x(:)));
    TS_DataMat = squashMe(TS_DataMat);
    for jo = 1:numGroups
        TS_DataMat(gi{jo},:) = squashMe(TS_DataMat(gi{jo},:)) + jo - 1;
    end
else
    numGroups = 0;
end

% --------------------------------------------------------------------------
%% Set the colormap
% --------------------------------------------------------------------------
if numGroups <= 1
    numColorMapGrads = 6; % number of gradations in each set of colourmap
    if strcmp(customColorMap,'redyellowblue')
        customColorMap = flipud(BF_GetColorMap('redyellowblue',numColorMapGrads,0));
    else
        customColorMap = gray(numColorMapGrads);
    end
elseif numGroups == 2
	% Special case to make a nice red and blue one
	customColorMap = [flipud(BF_GetColorMap('blues',9,0));flipud(BF_GetColorMap('reds',9,0))];
else
	% Use the same colors as GiveMeColors, but add brightness gradations to indicate magnitude
    numColorMapGrads = 20; % number of brightness gradations in each set of colourmap
    colormapBase = GiveMeColors(numGroups);
    customColorMap = [];
    for i = 1:numGroups
        customColorMap = [customColorMap; BF_MakeBrightenedColorMap(colormapBase{i},numColorMapGrads)];
    end
end

%-------------------------------------------------------------------------------
% Plotting
%-------------------------------------------------------------------------------
F = figure('color','w');

% ------------------------------------------------------------------------------
%% Plot the data matrix
% ------------------------------------------------------------------------------
colormap(customColorMap)
load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/stage_ordered(439)')

% TS_DataMat = TS_DataMat(stage_ordered,:);
load('HCTSA_N.mat', 'op_clust')
TS_DataMat = TS_DataMat(:,op_clust.ord)';
imagesc(TS_DataMat);

% ------------------------------------------------------------------------------
% Superimpose colored rectangles over NaN values
% ------------------------------------------------------------------------------
if colorNaNs && any(isnan(TS_DataMat(:)))
    [theNaNs_i,theNaNs_j] = find(isnan(TS_DataMat));
    fprintf(1,['Superimposing black rectangles over all %u NaNs in ' ...
                                'the data matrix\n'],length(theNaNs_i));
    for i = 1:length(theNaNs_i)
        rectangle('Position',[theNaNs_j(i)-0.5,theNaNs_i(i)-0.5,1,1],'FaceColor','k', ...
                        'EdgeColor','k')
    end
end

% --------------------------------------------------------------------------
%% Format the axes
% --------------------------------------------------------------------------
% Axis labels:
ax2 = gca;
ax2.FontSize = 8; % small font size (often large datasets)
ax2.TickLabelInterpreter = 'none'; % prevent displaying underscores as subscripts

% Rows: time series
% ax2.YTick = 1:numTS;
% ax2.YLim = [0.5,numTS+0.5];

% Columns: operations:
xlabel('Time (hours)')
% ax2.XLim = [0.5,numOps+0.5];
if numOps < 1000 % if too many operations, it's too much to list them all...
    ax2.XTick = 1:numOps;
    ax2.XTickLabel = Operations.Name;
    ax2.XTickLabelRotation = 90;
end

label_p = ylabel('Features','fontsize',17);
label_p.Position(2) = 3000; 
label_p.Position(1) = -50; 

% Add a color bar:

% cB = colorbar('southoutside');
% cB.Position = [0.130357142857143,0.188652866680432,0.183950910309393,0.029626436781609];  % Change last digit for the height of colorbar
% cB.Label.String = 'Normalised feature value';


cB.Ticks = 0:0.2:1;
cB.TickLabelInterpreter = 'none';


% Mark channel conditions
% yline(1,'k-','EEG','LineWidth',2,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',15,'FontWeight','bold')
% yline(6006,'k-','EOG','LineWidth',2,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',15,'FontWeight','bold')
% yline(12012,'k-','EMG','LineWidth',2,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',15,'FontWeight','bold')


end
