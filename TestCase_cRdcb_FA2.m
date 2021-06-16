fclose all;
close all;
clear;
clc;

%% Experimental
% Load data
load('Data/cRDCB FA2 Trac-Sep/cRdcb_FA2_ExpTracSep_Corridors.mat')
load('Data/cRDCB FA2 Trac-Sep/cRdcb_FA2_ExpTracSep.mat')

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
xlim([0,0.01])
ylim([0,45])
xlabel('Separation (mm)')
ylabel('Traction (N/mm^2)')

title('cRDCB FA2 Experimental Traction-Separation')

%% Trapezoidal Demo
load('Data/cRDCB FA2 Trac-Sep\cRdcb_FA2_TrapezoidTracSep.mat')

% Run ARCGen using ellipsoids
[charAvg, innerCorr, outerCorr,quadrants] = ARCGen_Rectangle(responseCurves,...
    'Diagnostics','on',...
    'nResamplePoints',75);