fclose all;
close all;
clear;
clc;

nResample = 100;
smfact = 1;

xlimits = [0,40];
ylimits = [0,120];

%% Generate Force-disp corridors without magnitude normalization
% Load data
load('Data/DCB Force-Disp/Dcb_ForceDisp.mat')

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
[~,~,segs] = selfintersect(outCorrNoNorm(:,1),outCorrNoNorm(:,2));
numbers = [1:1:length(outCorrNoNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
outCorrNoNorm = outCorrNoNorm(~any(indices,2),:);
clear indices

figure('Name','Force-Disp - No Normalization');
hold on;
for iPlot = 1:length(proCurveDataNoNorm)
    pExp = plot(proCurveDataNoNorm(iPlot).data(:,1),...
        proCurveDataNoNorm(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end


pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate Force-disp corridors with magnitude normalization
% Load data
load('Data/DCB Force-Disp/Dcb_ForceDisp.mat')

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
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

[~,~,segs] = selfintersect(outCorrNorm(:,1),outCorrNorm(:,2));
numbers = [1:1:length(outCorrNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
outCorrNorm = outCorrNorm(~any(indices,2),:);
clear indices

figure('Name','Force-Disp - Normalization');
hold on;
for iPlot = 1:length(proCurveDataNorm)
    pExp = plot(proCurveDataNorm(iPlot).data(:,1),...
        proCurveDataNorm(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate Force-disp corridors Remove Extaneous
% Load data
load('Data/DCB Force-Disp/Dcb_ForceDisp.mat')

[charAvgRemove, innCorrRemove, outCorrRemove,proCurveDataRemove] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'RemoveExtraneous',...
    'DeviationFact',   2);

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrRemove(:,1),innCorrRemove(:,2));
numbers = [1:1:length(innCorrRemove)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrRemove = innCorrRemove(~any(indices,2),:);
clear indices
[~,~,segs] = selfintersect(outCorrRemove(:,1),outCorrRemove(:,2));
numbers = [1:1:length(outCorrRemove)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
outCorrRemove = outCorrRemove(~any(indices,2),:);
clear indices

figure('Name','Force-Disp - Removed Extraneous');
hold on;
for iPlot = 1:length(proCurveDataRemove)
    pExp = plot(proCurveDataRemove(iPlot).data(:,1),...
        proCurveDataRemove(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvgRemove(:,1),charAvgRemove(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrRemove(:,1),innCorrRemove(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrRemove(:,1),outCorrRemove(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate Force-disp corridors crop to deviation factor
% Load data
load('Data/DCB Force-Disp/Dcb_ForceDisp.mat')

[charAvgLongCrop, innCorrLongCrop, outCorrLongCrop,proCurveDataLongCrop] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'CropToDeviationFactor',...
    'DeviationFact',   2);

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrLongCrop(:,1),innCorrLongCrop(:,2));
numbers = [1:1:length(innCorrLongCrop)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrLongCrop = innCorrLongCrop(~any(indices,2),:);
clear indices
[~,~,segs] = selfintersect(outCorrLongCrop(:,1),outCorrLongCrop(:,2));
numbers = [1:1:length(outCorrLongCrop)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
outCorrLongCrop = outCorrLongCrop(~any(indices,2),:);
clear indices

figure('Name','Force-Disp - Deviation Factor Crop');
hold on;
for iPlot = 1:length(proCurveDataLongCrop)
    pExp = plot(proCurveDataLongCrop(iPlot).data(:,1),...
        proCurveDataLongCrop(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvgLongCrop(:,1),charAvgLongCrop(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrLongCrop(:,1),innCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrLongCrop(:,1),outCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate R-Curve corridors without magnitude normalization
xlimits = [0,75];
ylimits = [0,1.2];
clear responseCurves
% Load data

load('Data/DCB Force-Disp/Dcb_rCurves.mat')

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off');

% Hackish way to clean up corridors by removing self-intersections
% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrNoNorm(:,1),innCorrNoNorm(:,2));
numbers = [1:1:length(innCorrNoNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrNoNorm = innCorrNoNorm(~any(indices,2),:);
clear indices
[~,~,segs] = selfintersect(outCorrNoNorm(:,1),outCorrNoNorm(:,2));
numbers = [1:1:length(outCorrNoNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
outCorrNoNorm = outCorrNoNorm(~any(indices,2),:);
clear indices


figure('Name','R-Curve - No Normalization');
hold on;
for iPlot = 1:length(proCurveDataNoNorm)
    pExp = plot(proCurveDataNoNorm(iPlot).data(:,1),...
        proCurveDataNoNorm(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end


pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate R-Curve corridors with magnitude normalization
% Load data
load('Data/DCB Force-Disp/Dcb_rCurves.mat')

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
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

[~,~,segs] = selfintersect(outCorrNorm(:,1),outCorrNorm(:,2));
numbers = [1:1:length(outCorrNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
outCorrNorm = outCorrNorm(~any(indices,2),:);
clear indices


figure('Name','R-Curve - Normalization');
hold on;
for iPlot = 1:length(proCurveDataNorm)
    pExp = plot(proCurveDataNorm(iPlot).data(:,1),...
        proCurveDataNorm(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate R-Curve corridors Remove Extaneous
% Load data
load('Data/DCB Force-Disp/Dcb_rCurves.mat')

[charAvgRemove, innCorrRemove, outCorrRemove,proCurveDataRemove] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'RemoveExtraneous',...
    'DeviationFact',   2);

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrRemove(:,1),innCorrRemove(:,2));
numbers = [1:1:length(innCorrRemove)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrRemove = innCorrRemove(~any(indices,2),:);
clear indices
[~,~,segs] = selfintersect(outCorrRemove(:,1),outCorrRemove(:,2));
numbers = [1:1:length(outCorrRemove)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
outCorrRemove = outCorrRemove(~any(indices,2),:);
clear indices

figure('Name','R-Curve - Removed Extraneous');
hold on;
for iPlot = 1:length(proCurveDataRemove)
    pExp = plot(proCurveDataRemove(iPlot).data(:,1),...
        proCurveDataRemove(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvgRemove(:,1),charAvgRemove(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrRemove(:,1),innCorrRemove(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrRemove(:,1),outCorrRemove(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate R-Curve corridors crop to deviation factor
% Load data
load('Data/DCB Force-Disp/Dcb_rCurves.mat')

[charAvgLongCrop, innCorrLongCrop, outCorrLongCrop,proCurveDataLongCrop] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'CropToDeviationFactor',...
    'DeviationFact',   0.75);

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrLongCrop(:,1),innCorrLongCrop(:,2));
numbers = [1:1:length(innCorrLongCrop)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrLongCrop = innCorrLongCrop(~any(indices,2),:);
clear indices

figure('Name','R-Curve - Deviation Factor Crop');
hold on;
for iPlot = 1:length(proCurveDataLongCrop)
    pExp = plot(proCurveDataLongCrop(iPlot).data(:,1),...
        proCurveDataLongCrop(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvgLongCrop(:,1),charAvgLongCrop(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrLongCrop(:,1),innCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrLongCrop(:,1),outCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate R-Curve corridors crop to deviation factor
% Load data
load('Data/DCB Force-Disp/Dcb_rCurves.mat')

[charAvgLongCrop, innCorrLongCrop, outCorrLongCrop,proCurveDataLongCrop] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics',     'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers',  'WeightedAverage',...
    'DeviationFact',   2);

% Hackish way to clean up corridors by removing self-intersections
% [~,~,segs] = selfintersect(innCorrLongCrop(:,1),innCorrLongCrop(:,2));
% numbers = [1:1:length(innCorrLongCrop)]';
% for iSeg = 1:size(segs,1)
%     indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
% end
% innCorrLongCrop = innCorrLongCrop(~any(indices,2),:);
% clear indices

figure('Name','R-Curve - Weighted Average');
hold on;
for iPlot = 1:length(proCurveDataLongCrop)
    pExp = plot(proCurveDataLongCrop(iPlot).data(:,1),...
        proCurveDataLongCrop(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvgLongCrop(:,1),charAvgLongCrop(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrLongCrop(:,1),innCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrLongCrop(:,1),outCorrLongCrop(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')