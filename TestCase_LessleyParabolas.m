fclose all;
close all;
clear;
clc;

%% Three parabolas used in paper
% Load data
load('Data/Lessley Parabolas/Lessley_Parabola_Processed.mat')
invalidCurves = [4];
validIndices = not([1:length(responseCurves)]==invalidCurves);
responseCurves = responseCurves(validIndices);

figure();
hold on;
cmap = cbrewer2('Paired',length(responseCurves))
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

load('Data/Lessley Parabolas/Lessley_ValidParabolas_Corridors.mat')
pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Characteristic Avg','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pInner = plot(innerCorr(:,1),innerCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
pOuter = plot(outerCorr(:,1),outerCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pInner])
xlim([0,1.6])
ylim([0,160])
title('Lessley et al. (2004) Parabolas')

%% Barabola with Outlier
% Load data
load('Data/Lessley Parabolas/Lessley_Parabola_Processed.mat')
load('Data/Lessley Parabolas/Lessley_OutlierParabola_Corridors.mat')

figure();
hold on;
cmap = cbrewer2('Paired',length(responseCurves))
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Characteristic Avg','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pInner = plot(innerCorr(:,1),innerCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
pOuter = plot(outerCorr(:,1),outerCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pInner])
xlim([0,2.2])
ylim([0,460])
title('Lessley et al. (2004) Parabolas with Outlier')