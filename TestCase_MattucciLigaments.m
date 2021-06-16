fclose all;
close all;
clear;
clc;

%% Mattucci Anterior Longitudinal Ligament
% Load data
load('Data/Mattucci Ligament Data/Mattucci_AnteriorLongitudinalLigament_QuasiStatic_NoFailure.mat')
load('Data/Mattucci Ligament Data/Mattucci_AnteriorLongitudinalLigament_QuasiStatic_NoFailure_Corridor.mat')

figure();
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','ARCGen - Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pArc = plot(innerCorr(:,1),innerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','ARCGen - Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outerCorr(:,1),outerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Outer - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

% Mattucci Parameters

C3 = 1839710;
C4 = 1.79e-05;
m = 1.167;
dtoe = 1.79;
xx = linspace(0,dtoe,20)';

mattucciAvg = [xx,C3.*(exp(C4.*xx)-1).*xx.^m];

fD = 116.0;
C5 = 140.4;
dLin = 2.7;
xx = linspace(dtoe,dLin,10)'

mattucciAvg = [mattucciAvg;[xx,fD+(xx-dtoe)*C5]];

A = 0.47;
B = -48.37;
C = 390.97;
D = -468.41;
dfail = 4.12;
xx = linspace(dLin,dfail,20)';

mattucciAvg = [mattucciAvg;[xx, A.*xx.^3 + B.*xx.^2 + C.*xx + D]];

pMattucci = plot(mattucciAvg(:,1),mattucciAvg(:,2),...
    'LineWidth',2.5,'Color',[55,126,184]./255,...
    'DisplayName','Char. Avg. - Mattucci');


legend([pExp,pAvg,pArc,pMattucci],'Location','best')
xlim([0,6.0])
ylim([0,700])

title('Mattucci - Anterior Longitudinal - Quasi-static')

%% Mattucci Posterior Longitudinal Ligament
% Load data
load('Data/Mattucci Ligament Data/Mattucci_PosteriorLongitudinalLigament_QuasiStatic_NoFailure.mat')
load('Data/Mattucci Ligament Data/Mattucci_PosteriorLongitudinalLigament_QuasiStatic_NoFailure_Corridor.mat')

figure();
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','ARCGen - Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pArc = plot(innerCorr(:,1),innerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','ARCGen - Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outerCorr(:,1),outerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Outer - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

% Mattucci Parameters

C3 = 3353314;
C4 = 3.49e-5;
m = 1.167;
dtoe = 0.79;
xx = linspace(0,dtoe,20)';

mattucciAvg = [xx,C3.*(exp(C4.*xx)-1).*xx.^m];

fD = 69;
C5 = 217.6;
dLin = 1.26;
xx = linspace(dtoe,dLin,10)';

mattucciAvg = [mattucciAvg;[xx,fD+(xx-dtoe)*C5]];

A = 0.63;
B = -72.94;
C = 397.73;
D = -215.29;
dfail = 2.69;
xx = linspace(dLin,dfail,20)';

mattucciAvg = [mattucciAvg;[xx, A.*xx.^3 + B.*xx.^2 + C.*xx + D]];

pMattucci = plot(mattucciAvg(:,1),mattucciAvg(:,2),...
    'LineWidth',2.5,'Color',[55,126,184]./255,...
    'DisplayName','Mattucci - Avg.');


legend([pExp,pAvg,pArc,pMattucci],'Location','best')
xlim([0,5.0])
ylim([0,750])

title('Mattucci - Posterior Longitudinal - Quasi-static')

%% Mattucci Ligamentum Flavum 
% Load data
load('Data/Mattucci Ligament Data/Mattucci_LigamentumFlavum_QuasiStatic_NoFailure.mat')
load('Data/Mattucci Ligament Data/Mattucci_LigamentumFlavum_QuasiStatic_NoFailure_Corridor.mat')

figure();
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','ARCGen - Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pArc = plot(innerCorr(:,1),innerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','ARCGen - Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outerCorr(:,1),outerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Outer - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

% Mattucci Parameters

C3 = 2513946;
C4 = 3.07e-7;
m = 2.761;
dtoe = 3.82;
xx = linspace(0,dtoe,20)';

mattucciAvg = [xx,C3.*(exp(C4.*xx)-1).*xx.^m];

fD = 119.7;
C5 = 117.7;
dLin = 4.29;
xx = linspace(dtoe,dLin,10)';

mattucciAvg = [mattucciAvg;[xx,fD+(xx-dtoe)*C5]];

A = -0.10;
B = -46.69;
C = 523.98;
D = -1206.13;
dfail = 5.53;
xx = linspace(dLin,dfail,20)';

mattucciAvg = [mattucciAvg;[xx, A.*xx.^3 + B.*xx.^2 + C.*xx + D]];

pMattucci = plot(mattucciAvg(:,1),mattucciAvg(:,2),...
    'LineWidth',2.5,'Color',[55,126,184]./255,...
    'DisplayName','Mattucci - Avg.');


legend([pExp,pAvg,pArc,pMattucci],'Location','best')
xlim([0,7.5])
ylim([0,600])

title('Mattucci - Ligamentum Flavum - Quasi-static')

%% Mattucci Capsular Ligament
% Load data
load('Data/Mattucci Ligament Data/Mattucci_CapsularLigament_QuasiStatic_NoFailure.mat')
load('Data/Mattucci Ligament Data/Mattucci_CapsularLigament_QuasiStatic_NoFailure_Corridor.mat')

figure();
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','ARCGen - Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pArc = plot(innerCorr(:,1),innerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','ARCGen - Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outerCorr(:,1),outerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Outer - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

% Mattucci Parameters

C3 = 1726474;
C4 = 1.55e-5;
m = 1.184;
dtoe = 1.33;
xx = linspace(0,dtoe,20)';

mattucciAvg = [xx,C3.*(exp(C4.*xx)-1).*xx.^m];

fD = 50.2;
C5 = 82.2;
dLin = 1.76;
xx = linspace(dtoe,dLin,10)';

mattucciAvg = [mattucciAvg;[xx,fD+(xx-dtoe)*C5]];

A = -0.17;
B = -14.72;
C = 135.58;
D = -106.90;
dfail = 4.16;
xx = linspace(dLin,dfail,20)';

mattucciAvg = [mattucciAvg;[xx, A.*xx.^3 + B.*xx.^2 + C.*xx + D]];

pMattucci = plot(mattucciAvg(:,1),mattucciAvg(:,2),...
    'LineWidth',2.5,'Color',[55,126,184]./255,...
    'DisplayName','Mattucci - Avg.');


legend([pExp,pAvg,pArc,pMattucci],'Location','best')
xlim([0,7.0])
ylim([0,350])

title('Mattucci - Capsular Ligament - Quasi-static')

%% Mattucci Interspinous Ligament
% Load data
load('Data/Mattucci Ligament Data/Mattucci_InterspinousLigament_QuasiStatic_NoFailure.mat')
load('Data/Mattucci Ligament Data/Mattucci_InterspinousLigament_QuasiStatic_NoFailure_Corridor.mat')

figure();
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','ARCGen - Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',0.0.*[1,1,1]);
pArc = plot(innerCorr(:,1),innerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','ARCGen - Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outerCorr(:,1),outerCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Outer - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

% Mattucci Parameters

C3 = 610863;
C4 = 4.86e-6;
m = 1.050;
dtoe = 3.02;
xx = linspace(0,dtoe,20)';

mattucciAvg = [xx,C3.*(exp(C4.*xx)-1).*xx.^m];

fD = 28.6;
C5 = 19.4;
dLin = 3.60;
xx = linspace(dtoe,dLin,10)';

mattucciAvg = [mattucciAvg;[xx,fD+(xx-dtoe)*C5]];

A = 0.32;
B = -8.06;
C = 64.90;
D = -104.41;
dfail = 6.93;
xx = linspace(dLin,dfail,20)';

mattucciAvg = [mattucciAvg;[xx, A.*xx.^3 + B.*xx.^2 + C.*xx + D]];

pMattucci = plot(mattucciAvg(:,1),mattucciAvg(:,2),...
    'LineWidth',2.5,'Color',[55,126,184]./255,...
    'DisplayName','Mattucci - Avg.');

legend([pExp,pAvg,pArc,pMattucci],'Location','best')
xlim([0,10])
ylim([0,150])

title('Mattucci - Interspinous Ligament - Quasi-static')