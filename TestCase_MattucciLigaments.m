fclose all;
% close all;
clear;
clc;

nResample = 150;
nCorrPts = 400;

xlimits = [0,6.0];
ylimits = [0,700];
%% Generate Anterior Longitudinal, No Normalization
% Load data
load('Data/Mattucci Ligament Data/Mattucci_AnteriorLongitudinalLigament_QuasiStatic_NoFailure')
invalidCurves = [11,14,15,16]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1);

figure('Name','Anterior Longitudinal- No Normalization');
title('Anterior Longitudinal')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
xx = linspace(dtoe,dLin,10)' ;
 
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

pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')

%% Generate Anterior Longitudinal, Normalized
% Load data
load('Data/Mattucci Ligament Data/Mattucci_AnteriorLongitudinalLigament_QuasiStatic_NoFailure')
invalidCurves = [11,14,15,16]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'on',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1,...
    'nWarpCtrlPts',     0,...
    'warpingPenalty',   0);

figure('Name','Anterior Longitudinal- Normalized');
title('Anterior Longitudinal')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
xx = linspace(dtoe,dLin,10)' ;
 
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

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

xlimits = [0,5.0];
ylimits = [0,600];
%% Generate Posterior Longitudinal, No Normalization
% Load data
load('Data/Mattucci Ligament Data/Mattucci_PosteriorLongitudinalLigament_QuasiStatic_NoFailure')
invalidCurves = [5,6,14]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1);

figure('Name','Posterior Longitudinal- No Normalization');
title('Posterior Longitudinal')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate Posterior Longitudinal, Normalized
% Load data
load('Data/Mattucci Ligament Data/Mattucci_PosteriorLongitudinalLigament_QuasiStatic_NoFailure')
invalidCurves = [5,6,14]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'on',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1,...
    'nWarpCtrlPts',     0,...
    'warpingPenalty',   1e-3);

figure('Name','Posterior Longitudinal- Normalized');
title('Posterior Longitudinal')
hold on;
for iPlot = 1:length(proCurveDataNorm)
    pExp = plot(proCurveDataNorm(iPlot).data(:,1),...
        proCurveDataNorm(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

xlimits = [0,7.5];
ylimits = [0,600];
%% Generate Ligamentum Flavum, No Normalization
% Load data
load('Data/Mattucci Ligament Data/Mattucci_LigamentumFlavum_QuasiStatic_NoFailure')

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1);

figure('Name','Ligamentum Flavum- No Normalization');
title('Ligamentum Flavum')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate Ligamentum Flavum, Normalized
% Load data
load('Data/Mattucci Ligament Data/Mattucci_LigamentumFlavum_QuasiStatic_NoFailure')

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1,...
    'nWarpCtrlPts',     0,...
    'warpingPenalty',   1e-3);

figure('Name','Ligamentum Flavum- Normalized');
title('Ligamentum Flavum')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

xlimits = [0,7.5];
ylimits = [0,400];
%% Generate Capsular Ligament, No Normalization
% Load data
load('Data/Mattucci Ligament Data/Mattucci_CapsularLigament_QuasiStatic_NoFailure')
invalidCurves = [1,6,9]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1);

figure('Name','Capsular Ligament- No Normalization');
title('Capsular Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate Capsular Ligament, Normalized
% Load data
load('Data/Mattucci Ligament Data/Mattucci_CapsularLigament_QuasiStatic_NoFailure')
invalidCurves = [1,6,9]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1,...
    'nWarpCtrlPts',     0,...
    'warpingPenalty',   1e-3);

figure('Name','Capsular Ligament- Normalized');
title('Capsular Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

xlimits = [0,12];
ylimits = [0,150];
%% Generate Interspinous Ligament, No Normalization
% Load data
load('Data/Mattucci Ligament Data/Mattucci_InterspinousLigament_QuasiStatic_NoFailure')
invalidCurves = [9,10,13]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNoNorm, innCorrNoNorm, outCorrNoNorm,proCurveDataNoNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'off',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1);

figure('Name','Interspinous Ligament- No Normalization');
title('Interspinous Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNoNorm(:,1),charAvgNoNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNoNorm(:,1),innCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNoNorm(:,1),outCorrNoNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')

%% Generate Interspinous Ligament, Normalized
% Load data
load('Data/Mattucci Ligament Data/Mattucci_InterspinousLigament_QuasiStatic_NoFailure')
invalidCurves = [9,10,13]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

[charAvgNorm, innCorrNorm, outCorrNorm,proCurveDataNorm] = ...
    ARCGen_Ellipse(responseCurves,...
    'Diagnostics', 'off',...
    'nResamplePoints', nResample,...
    'NormalizeCurves', 'on',...
    'handleOutliers', 'off',...
    'CorridorRes',      nCorrPts,....
    'EllipseKFact',     1,...
    'nWarpCtrlPts',     0,...
    'warpingPenalty',   0);

figure('Name','Interspinous Ligament- Normalized');
title('Interspinous Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

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
    'DisplayName','Char. Avg. - Mattucci'); 

pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim(xlimits)
ylim(ylimits)
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')