
%%%%%%% Had to change l.114: Change TimeSeries.Data into cell array 

function TS_PlotDataMatrix_edited(varargin)
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

%%%% Only EEG (439)

numTS = 1166;
TimeSeries = TimeSeries(1:numTS,:);
TS_DataMat = TS_DataMat(1:numTS,:);

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
f = figure('color','w');

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

label_p = ylabel('Operations','fontsize',17);
label_p.Position(2) = 3000; 
label_p.Position(1) = -45; 

% Add a color bar:
cB = colorbar('eastoutside');
cB.Position = [0.957165258576549,0.4701532324569,0.007680491551459,0.360683760683761];  % Change last digit for the height of colorbar
cB.Label.String = 'Feature value';


cB.Ticks = 0:0.2:1;
cB.TickLabelInterpreter = 'none';


% Mark channel conditions
% yline(1,'k-','EEG','LineWidth',2,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',15,'FontWeight','bold')
% yline(6006,'k-','EOG','LineWidth',2,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',15,'FontWeight','bold')
% yline(12012,'k-','EMG','LineWidth',2,'LabelHorizontalAlignment','left','LabelVerticalAlignment','bottom','FontSize',15,'FontWeight','bold')

pos = get(f,'position');

% %% To inspect a couple of TS
% figure;  
% [ha, pos] = tight_subplot(2,2,[.090 .1],[.1 .05],[.05 .05]);
% ax = gca;
% 
% corr = 230;
% fal = 573;
% 
% axes(ha(1)); plot(TimeSeries.Data{corr,:})
% axes(ha(2)); plot(TimeSeries.Data{fal,:})
% axes(ha(3)); imagesc(TS_DataMat(:,corr));
% axes(ha(4)); imagesc(TS_DataMat(:,fal));
% numColorMapGrads = 6; 
% customColorMap = flipud(BF_GetColorMap('redyellowblue',numColorMapGrads,0));
% colormap(customColorMap)

%%%%%%%%%%%%%%%%%
%%%% Agreed epochs
%%%%%%%%%%%%%%%%%

%%%% Same signature
SS_wake = 86;
SS_wake_2 = 1040;
SS_N2 = 476;
SS_N2_2 = 162;
SS_N3 = 195;
SS_N3_2 = 710;

%%%% No signature
NS_N1 = 274;
NS_N1_2 = 275;
NS_N2 = 285;
NS_N2_2 = 298;
NS_R = 1119;
NS_R_2 = 624;


%%%%%%%%%%%%%%%%%
%%%% Misclassified epochs
%%%%%%%%%%%%%%%%%

%%%% Same signature
SS_WN1 = 1049;
SS_WR= 1051;
SS_N2R = 670;
SS_N2N3 = 697;
SS_N3R = 761;
SS_N3W = 405;


%%%% No signature
NS_N1W = 157;
NS_N1N2 = 647;
NS_N2R = 461;
NS_N2N3 = 327 ;
NS_RW = 945;
NS_RN1 = 918;


% Epochs_SS = [SS_wake SS_WN1 SS_wake_2 SS_WR SS_N2 SS_N2R SS_N2_2 SS_N2N3 SS_N3 SS_N3R SS_N3_2 SS_N3W ];
Epochs_NS = [SS_wake SS_WN1 SS_wake_2 SS_WR  NS_N1 NS_N1W NS_N1_2 NS_N1N2 NS_N2 NS_N2R NS_N2_2 NS_N2N3 NS_R NS_RW NS_R_2 NS_RN1];
 
TITLE = {'wake','N1','wake','REM','N2','REM','N2','N3','N3','REM','N3','wake'};
TITLE_NS = {'wake','N1','wake','REM','N1','wake','N1','N2','N2','REM','N2','N3','REM','wake','REM','N1'};
Epochs_NS_str = {'SS_wake' 'SS_WN1' 'SS_WR' 'SS_wake_2'  'NS_N1' 'NS_N1W' 'NS_N1_2' 'NS_N1N2' 'NS_N2' 'NS_N2R' 'NS_N2_2' 'NS_N2N3' 'NS_R' 'NS_RW' 'NS_R_2' 'NS_RN1'};


% load('HCTSA_N.mat', 'TimeSeries')

%%%%%% Plot TS
for subplot = 1:16

    TS = TimeSeries.Data{Epochs_NS(subplot),:};
    
    f = figure; ax=gca; 
    plot(TS,'LineWidth',1.5)
    % plot(TS)

    title(TITLE_NS{subplot},'FontSize',55)
    xticks([0:1280:3840]);
    ax.XTickLabels = {'0','','','30'};
    
    % Label on the right axis for the  plotson the right
    if subplot == 2 || subplot ==4  || subplot ==6 || subplot == 8 || subplot ==10 || subplot ==12 || subplot ==14 || subplot ==16
        ax.YAxisLocation = 'right';
        ax.Position = [0.04,0.23,0.85,0.55];
    else
        ax.Position = [0.090,0.23,0.85,0.55];
    end
    
    % Y axis limit for wake
    if subplot == 1 || subplot == 2  || subplot == 3 || subplot == 4
         ylim([-0.15 0.15])
         yticks(-0.15:0.15:0.15);
    % Y axis limit for N1
    elseif subplot == 5 || subplot == 6 || subplot == 7 || subplot == 8
        ylim([-0.40 0.2])
        yticks(-0.40:0.2:0.2);
    % Y axis limit for N2
    elseif subplot == 9 || subplot == 10 || subplot == 11 || subplot == 12
         ylim([-0.15 0.15])
         yticks(-0.15:0.15:0.15);
    % Y axis limit for REM
    elseif subplot == 13 || subplot == 14 || subplot == 15 || subplot == 16
         ylim([-0.15 0.15])
         yticks(-0.15:0.15:0.15);
    end
    
    f.Position = [1,377,1449,275];
    set(gca,'box','off') 
    ax.FontSize = 45;

    set(f, 'Color', 'w')
    fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms/EachTS';
    export_fig([fpath filesep sprintf('%s',Epochs_NS_str{subplot})],'-r 300')

end 

%%%%%% Plot value colorbars
for subplot = 1:16
    
    g = figure; ax = gca; 
    imagesc(TS_DataMat(:,Epochs_NS(subplot)));

    numColorMapGrads = 6; 
    customColorMap = flipud(BF_GetColorMap('redyellowblue',numColorMapGrads,0));
    colormap(customColorMap)
    
    set(gca,'xtick',[]); set(gca,'ytick',[])
    set(gca,'visible','off')
    g.Position = [274,332,97,188];
    set(g, 'Color', 'w')
    
    fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms/EachFV';
    export_fig([fpath filesep sprintf('%s',Epochs_NS_str{subplot})],'-r 300')

end

set(a, 'Color', 'w')
% fpath = '/Users/nico/Documents/HCTSA/Analysis/hypnograms';
% export_fig([fpath filesep 'hypno_EEG_cl(439)_values_NS'],'-r 300')


disp('ok')


%% Make the selection

load('/Users/nico/Documents/HCTSA/Analysis/hypnograms/statsOut_allepochs_3ch(439)')
original_labels = statsOut.scoredTest;
cluster_decision = statsOut.predictTest;

% Give equivalent label for original labels
a = find(original_labels == 0);
b = find(original_labels == 1);
c = find(original_labels == 2);
d = find(original_labels == 3);
e = find(original_labels == 5);

% Give equivalent label for original labels
f = find(cluster_decision == 0);
g = find(cluster_decision == 1);
h = find(cluster_decision == 2);
i = find(cluster_decision == 3);
j = find(cluster_decision == 5);


False = a(find(cluster_decision(a) ~= 0));
false_ST = cluster_decision(False);


% load('HCTSA_N.mat','TimeSeries')

STG = {'Wake','N1','N2','N3','','REM'};

%for subplot = 1:length(False)
for subplot = 80:150

    f = figure; ax=gca;
    TS = TimeSeries.Data{False(subplot),:};

    plot(TS)
    
    xticks([0:1280:3840]);
    xticklabels({'0','10','20','30'})
    ylim([-0.15 0.15])
    yticks(-0.15:0.15:0.15);


    The_stage = cellstr(STG(false_ST(subplot)+1));
    The_ID = string(False(subplot));
    title([The_stage The_ID])
    
    ax.Position = [0.08,0.23,0.88,0.40];
    f.Position = [1,400,1449,275];
    set(gca,'box','off') 
    ax.FontSize = 45;
    
end 


%% CTL epoch

CTL = 139;

f = figure; ax=gca;
TS = TimeSeries.Data{CTL,:};

plot(TS)

xticks([0:1280:3840]);
xticklabels({'0','','','30'})
ylim([-0.15 0.15])
yticks(-0.15:0.15:0.15);

ax.Position = [0.08,0.23,0.88,0.40];
f.Position = [1,10,1449,275];
set(gca,'box','off') 
ax.FontSize = 45;

end




 

