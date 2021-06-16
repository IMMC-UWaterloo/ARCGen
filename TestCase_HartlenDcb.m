fclose all;
close all;
clear;
clc;

%% DCB force-displacement
% Load data
load('Data/DCB Force-Disp/Dcb_ForceDisp.mat')
load('Data/DCB Force-Disp/Dcb_ForceDisp_corridor.mat')

figure();
hold on;
cmap = cbrewer2('Paired',length(responseCurves));
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end


pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pCorr = plot(innerCorr(:,1),innerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outerCorr(:,1),outerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr],'Location','best')
xlim([0,40])
ylim([0,120])
xlabel('Displacement (mm)')
ylabel('Force (N)')

title('DCB Force-Displacement')

%% DCB rCurves
% Load data
load('Data/DCB Force-Disp/Dcb_rCurves.mat')
load('Data/DCB Force-Disp/Dcb_rCurves_corridor.mat')

figure();
hold on;
cmap = cbrewer2('Paired',length(responseCurves));
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end


pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pCorr = plot(innerCorr(:,1),innerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outerCorr(:,1),outerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr],'Location','best')
xlim([0,80])
ylim([0.2,1.1])
xlabel('Displacement (mm)')
ylabel('Force (N)')

title('DCBR-Curves')