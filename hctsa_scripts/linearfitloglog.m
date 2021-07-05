

%% Compute the linear fit to the log log (SP_Summaries_...) for each stage

% 1) Get the TS from each stage from 439 (average?)
% 2) Get the lf_a2 for each stage using SP_Summaries
% 3) Plot the log log spectrum with the linear fit (5 plots)


%%%%%% Get the TS for each stage (439)

% Load data (normalised FV)
sub = '439';
cd(sprintf('/Users/nico/Documents/MATLAB/hctsa-master/HCTSA_%s',sub))
load('HCTSA.mat','TimeSeries','TS_DataMat')

% Remove wakefulness before sleep
Sleep_inter = 151:1164;
TimeSeries = TimeSeries(Sleep_inter,:);
TS_DataMat = TS_DataMat(Sleep_inter,:);

% Get TS for each stage
wake = find(TimeSeries.Group == '0');
N1 = find(TimeSeries.Group == '1');
N2 = find(TimeSeries.Group == '2');
N3 = find(TimeSeries.Group == '3');
rem = find(TimeSeries.Group == '5');

stages = [{wake} {N1} {N2} {N3} {rem}];

% TS = time series for each stage
% mean_TS = average of data points for each stage
for S = 1:5
    TS{S,1} = cell2mat(TimeSeries.Data(stages{S})')';
    mean_TS(S,:) = mean(TS{S,1});
end


%%% Find representative TS for each stage (closest) to lf_a2

for S = 1:5
    
    for i = 1:size(TS{S},1)
        
        [out,~,~,~] = SP_Summaries_edited(TS{S,1}(i,:),'welch','rect');
        
        R_lf{S,1}(i,1) = out.linfitloglog_lf_a2;

    end
    
end

Aver_lf_a2 = [-1.5705 -1.6766 -2.0310 -2.8343 -1.8349];  % Average raw values

% To get representative TS, we find the TS whose lf_a2 is closest to
% average lf_a2
for S = 1:5
    [~,index] = min(abs(R_lf{S}-Aver_lf_a2(S)));
    REP_TS(S,:) = TS{S}(index,:);
end


%%%%%%%% Get the lf_a2 for each stage using SP_Summaries, using
%%%%%%%% representative TS
%%
for S = 1:5

    y = REP_TS(S,:);   % input TS
    [out,w,SS,r_lf] = SP_Summaries_edited(y,'welch','rect');
    
    lf_a2(S,1) = out.linfitloglog_lf_a2;
    power{S,1} = w;
    Freq{S,1} = SS;
    lowfreq{S,1} = r_lf;
end

%%%%%%% Plot log-log power spectrum for each stage
STAGE = {'Wake','N1','N2','N3','REM'};
%%
figure;ax=gca;

for S = 1:5
    
    w = power{S,1};
    SS = Freq{S,1};
    r_lf = lowfreq{S,1};
    
    subplot(1,5,S);

    plot(log(w(r_lf)),log(SS(r_lf)),'LineWidth',1);
    
    hold on
    x = -4:0.25:0.5; 

    m  = Aver_lf_a2(S);  % Specify slope
    x1 = -2.5; % Specify starting x
    y1 = -7;  % Specify starting y

    Y = m*(x - x1) + y1;
    hPlot = plot(x,Y,'LineWidth',2);
    
    title(sprintf('%s',STAGE{S}))
    
    ylim([-20 -2])
    xlim([-7 1])
    
    box off
    set(gcf,'color','white')


end

hold off

f.Position = [53,279,1388,518];

%%
%%%%%%%%% Plot the slope only

g = figure; 

for S = 1:5

    x = min(log(w(r_lf))):0.25:max(log(w(r_lf))); % Defines the domain as [-15,25] with a refinement of 0.25

    m  = lf_a2(S);  % Specify slope
    x1 = -3.2; % Specify starting x
    y1 = -7;  % Specify starting y

    Y = m*(x - x1) + y1;
    hPlot = plot(x,Y);
    
    title(sprintf('%s',STAGE{S}))
    
    g.Position = [462,313,406,420];
    
    hold on
    
end

legend('Wake','N1','N2','N3','REM')
