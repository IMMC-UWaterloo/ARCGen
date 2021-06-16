fclose all;
close all;
clear;
clc;

%% Three parabolas used in paper
% Load data
load('Data/Watson 7333 12R RDCB/WatsonRdcb_7333_12Thou.mat')
load('Data/Watson 7333 12R RDCB/WatsonRdcb_7333_12Thou_Corridors.mat')

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
xlim([0,0.4])
ylim([0,850])
xlabel('Displacement (mm)')
ylabel('Force (N)')

title('Watson et al (2018) RDCB - 7333')