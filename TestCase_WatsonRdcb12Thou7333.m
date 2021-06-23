fclose all;
close all;
clear;
clc;

nResample = 100;
smfact = 5;

xlimits = [0,0.40];
ylimits = [0,900];

%% Generate Force-Disp corridors without magnitude normalization
% Load data
load('Data/Watson 7333 12R RDCB/WatsonRdcb_7333_12Thou.mat')

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off');

figure('Name','Force-Disp - No Normalization');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Hackish way to clean up corridors by removing self-intersections
[~,~,segs] = selfintersect(innCorrNoNorm(:,1),innCorrNoNorm(:,2));
numbers = [1:1:length(innCorrNoNorm)]';
for iSeg = 1:size(segs,1)
    indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
end
innCorrNoNorm = innCorrNoNorm(~any(indices,2),:);
clear indices
% [~,~,segs] = selfintersect(outCorrNoNorm(:,1),outCorrNoNorm(:,2));
% numbers = [1:1:length(outCorrNoNorm)]';
% for iSeg = 1:size(segs,1)
%     indices(:,iSeg) = and(numbers>segs(iSeg,1), numbers<segs(iSeg,2));
% end
% outCorrNoNorm = outCorrNoNorm(~any(indices,2),:);
% clear indices

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

%% Generate Force-Disp corridors without magnitude normalization
% Load data
load('Data/Watson 7333 12R RDCB/WatsonRdcb_7333_12Thou.mat')

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
    ARCGen_Rectangle(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers', 'off');

figure('Name','Force-Disp - Normalization');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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