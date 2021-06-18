fclose all;
% close all;
clear;
clc;

nResample = 100;
smfact = 1;

xlimits = [0,4.5];
ylimits = [0,1200];

lobdellCorrs = readmatrix('Data/Kroell 1971 Thorax Response/Lobdell 16mph Corridors.csv');
charAvgLobdell = lobdellCorrs(:,1:2);
charAvgLobdell = [charAvgLobdell(:,1)+0.5, charAvgLobdell(:,2)-150];
innCorrLobdell = lobdellCorrs(~isnan(lobdellCorrs(:,3)),3:4);
innCorrLobdell = [innCorrLobdell(:,1)+0.5, innCorrLobdell(:,2)-150];
outCorrLobdell = lobdellCorrs(~isnan(lobdellCorrs(:,5)),5:6);
outCorrLobdell = [outCorrLobdell(:,1)+0.5, outCorrLobdell(:,2)-150];

%% Generate corridors without magnitude normalization
% Load data
load('Data/Kroell 1971 Thorax Response/KroellThoraxResponse_1971.mat')

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off');

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrNoNorm(:,1),innCorrNoNorm(:,2));
numbers = [1:1:length(innCorrNoNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrNoNorm = innCorrNoNorm(~any(indices,2),:);
clear indices


figure('Name','No Normalization');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvgLob = plot(charAvgLobdell(:,1),charAvgLobdell(:,2),'-',...
    'DisplayName','Char. Avg. - Lobdell','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrLob = plot(innCorrLobdell(:,1),innCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Corridors - Lobdell',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(outCorrLobdell(:,1),outCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Outer',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr,pAvgLob,pCorrLob], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate corridors with magnitude normalization
% Load data
load('Data/Kroell 1971 Thorax Response/KroellThoraxResponse_1971.mat')

[charAvgNorm, innCorrNorm, outCorrNorm ,proCurveDataNorm] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers', 'off');

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrNorm(:,1),innCorrNorm(:,2));
numbers = [1:1:length(innCorrNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrNorm = innCorrNorm(~any(indices,2),:);
clear indices


figure('Name','Normalization');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvgLob = plot(charAvgLobdell(:,1),charAvgLobdell(:,2),'-',...
    'DisplayName','Char. Avg. - Lobdell','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrLob = plot(innCorrLobdell(:,1),innCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Corridors - Lobdell',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(outCorrLobdell(:,1),outCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Outer',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr,pAvgLob,pCorrLob], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate corridors, removing extraneous curves
% Load data
load('Data/Kroell 1971 Thorax Response/KroellThoraxResponse_1971.mat')

[charAvgRemoved, innCorrRemoved, outCorrRemoved ,proCurveDataRemoved] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'RemoveExtraneous',...
    'DeviationFact',   2);

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrRemoved(:,1),innCorrRemoved(:,2));
numbers = [1:1:length(innCorrRemoved)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrRemoved = innCorrRemoved(~any(indices,2),:);
clear indices

figure('Name','Extraneous Curves Removed');
hold on;
for iPlot = 1:length(proCurveDataRemoved)
    pExp = plot(proCurveDataRemoved(iPlot).data(:,1),...
        proCurveDataRemoved(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvgLob = plot(charAvgLobdell(:,1),charAvgLobdell(:,2),'-',...
    'DisplayName','Char. Avg. - Lobdell','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrLob = plot(innCorrLobdell(:,1),innCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Corridors - Lobdell',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(outCorrLobdell(:,1),outCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Outer',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

pAvg = plot(charAvgRemoved(:,1),charAvgRemoved(:,2),'.-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrRemoved(:,1),innCorrRemoved(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrRemoved(:,1),outCorrRemoved(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr,pAvgLob,pCorrLob], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate corridors, cropping all curves to shortest
% Load data
load('Data/Kroell 1971 Thorax Response/KroellThoraxResponse_1971.mat')

[charAvgShortCrop, innCorrShortCrop, outCorrShortCrop ,proCurveDataShortCrop] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'CropToShortest',...
    'DeviationFact',   2);

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrShortCrop(:,1),innCorrShortCrop(:,2));
numbers = [1:1:length(innCorrShortCrop)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrShortCrop = innCorrShortCrop(~any(indices,2),:);
clear indices


figure('Name','Crop To Shortest');
hold on;
for iPlot = 1:length(proCurveDataShortCrop)
    pExp = plot(proCurveDataShortCrop(iPlot).data(:,1),...
        proCurveDataShortCrop(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvgLob = plot(charAvgLobdell(:,1),charAvgLobdell(:,2),'-',...
    'DisplayName','Char. Avg. - Lobdell','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrLob = plot(innCorrLobdell(:,1),innCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Corridors - Lobdell',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(outCorrLobdell(:,1),outCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Outer',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

pAvg = plot(charAvgShortCrop(:,1),charAvgShortCrop(:,2),'.-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrShortCrop(:,1),innCorrShortCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrShortCrop(:,1),outCorrShortCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr,pAvgLob,pCorrLob], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate corridors, cropping to a factor of median
% Load data
load('Data/Kroell 1971 Thorax Response/KroellThoraxResponse_1971.mat')

[charAvgLongCrop, innCorrLongCrop, outCorrLongCrop ,proCurveDataLongCrop] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'CropToDeviationFactor',...
    'DeviationFact',   1);

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrLongCrop(:,1),innCorrLongCrop(:,2));
numbers = [1:1:length(innCorrLongCrop)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrLongCrop = innCorrLongCrop(~any(indices,2),:);
clear indices


figure('Name','Crop To Deviation Factor');
hold on;
for iPlot = 1:length(proCurveDataLongCrop)
    pExp = plot(proCurveDataLongCrop(iPlot).data(:,1),...
        proCurveDataLongCrop(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvgLob = plot(charAvgLobdell(:,1),charAvgLobdell(:,2),'-',...
    'DisplayName','Char. Avg. - Lobdell','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrLob = plot(innCorrLobdell(:,1),innCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Corridors - Lobdell',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(outCorrLobdell(:,1),outCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Outer',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

pAvg = plot(charAvgLongCrop(:,1),charAvgLongCrop(:,2),'.-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrLongCrop(:,1),innCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrLongCrop(:,1),outCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr,pAvgLob,pCorrLob], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate corridors, Weighted Average
% Load data
load('Data/Kroell 1971 Thorax Response/KroellThoraxResponse_1971.mat')

[charAvgWeight, innCorrWeight, outCorrWeight ,proCurveDataWeight] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'WeightedAverage',...
    'DeviationFact',   3);

% Hackish way to clean up corridors by removing self-intersections of INNER
% ONLY
[~,~,segs] = selfintersect(innCorrWeight(:,1),innCorrWeight(:,2));
numbers = [1:1:length(innCorrWeight)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrWeight = innCorrWeight(~any(indices,2),:);


figure('Name','Weighted Averages');
hold on;
for iPlot = 1:length(proCurveDataWeight)
    pExp = plot(proCurveDataWeight(iPlot).data(:,1),...
        proCurveDataWeight(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvgLob = plot(charAvgLobdell(:,1),charAvgLobdell(:,2),'-',...
    'DisplayName','Char. Avg. - Lobdell','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrLob = plot(innCorrLobdell(:,1),innCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Corridors - Lobdell',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(outCorrLobdell(:,1),outCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Outer',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

pAvg = plot(charAvgWeight(:,1),charAvgWeight(:,2),'.-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrWeight(:,1),innCorrWeight(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrWeight(:,1),outCorrWeight(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr,pAvgLob,pCorrLob], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')
