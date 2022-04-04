%% NBDL Head Kinematics, 15g Frontal Decelerations, Head Displacements
% This script was used to generate plots for Hartlen & Cronin (2022). This
% script also serves as a great example of how the number of warping points
% and warping penalty factor influences the resulting corridors. The
% parameters provided with this script were found to produce the best
% results. 
%
% This script produces corridors for head displacements. Please refer to 
% "TestCase_NBDL_15gFrontal_Disp.m" for head accelerations. 

% Please refer to test case readme document for more details. 
%
% Dataset Citation:
%    Ewing, C. L., & Thomas, D. J. (1972). "Human Head and Neck Response to
%       Impact Acceleration." Naval Aerospace Medical Research Lab
%       Pensacola Fl.
%
%    National Highway Traffic Safety Administration. (2017). "Biomechanics
%       Test Database." 
%       https://www.nhtsa.gov/research-data/databases-and-software
%
% Copyright (c) 2022 Devon C. Hartlen 

%% MATLAB initialization
fclose all;
close all;
clear;
clc;

addpath('../'); % Ensure ARCGen is on execution path

% Set resample points and corridor resolution. 500 points will demonstrate
% runtime difference between `UseParallel` set to `on` and `off` without
% taking excessively long to run. 
nResample = 500;
nCorrPts = 500;

%% Generate Anterior Longitudinal Ligament Corridors
% Load data
load('Mattucci Ligament Data/Mattucci_AnteriorLongitudinalLigament_QuasiStatic_NoFailure')
% Exclude curves to match published Mattucci averaging technique
invalidCurves = [11,14,15,16]'; 
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

% Execute ARCGen. Do not output processed data structure. 
[charAvg, innCorr, outCorr, ~] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts);

figure('Name','Anterior Longitudinal- Normalized');
title('Anterior Longitudinal Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Construct Mattucci Average. Parameters from Mattucci & Cronin (2015)
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
 
% Plot Mattucci average
pMattucci = plot(mattucciAvg(:,1),mattucciAvg(:,2),... 
    'LineWidth',2.5,'Color',[55,126,184]./255,... 
    'DisplayName','Char. Avg. - Mattucci'); 

% Plot ARCGen average and corridors
pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim([0,6.0])
ylim([0,700])
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')


%% Generate Posterior Longitudinal Ligament Corridors
% Load data
load('Mattucci Ligament Data/Mattucci_PosteriorLongitudinalLigament_QuasiStatic_NoFailure')
invalidCurves = [5,6,14]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

% Execute ARCGen. Do not output processed data structure. 
[charAvg, innCorr, outCorr, ~] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts);

figure('Name','Posterior Longitudinal Ligament');
title('Posterior Longitudinal Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Construct Mattucci Average. Parameters from Mattucci & Cronin (2015)
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

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim([0,5.0])
ylim([0,600])
grid on
xlabel('Deflection (m)')
ylabel('Force (N)')


%% Generate Ligamentum Flavum Corridors
% Load data
load('Mattucci Ligament Data/Mattucci_LigamentumFlavum_QuasiStatic_NoFailure')

% Execute ARCGen. Do not output processed data structure. 
[charAvg, innCorr, outCorr, ~] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts);

figure('Name','Ligamentum Flavum');
title('Ligamentum Flavum')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Construct Mattucci Average. Parameters from Mattucci & Cronin (2015)
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

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim([0,7.5])
ylim([0,600])
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')


%% Generate Capsular Ligament Corridors
% Load data
load('Mattucci Ligament Data/Mattucci_CapsularLigament_QuasiStatic_NoFailure')
invalidCurves = [1,6,9]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

% Execute ARCGen. Do not output processed data structure. 
[charAvg, innCorr, outCorr, ~] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts);

figure('Name','Capsular Ligament');
title('Capsular Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Construct Mattucci Average. Parameters from Mattucci & Cronin (2015)
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

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim([0,7.5])
ylim([0,400])
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')


%% Generate Interspinous Ligament Corridors
% Load data
load('Mattucci Ligament Data/Mattucci_InterspinousLigament_QuasiStatic_NoFailure')
invalidCurves = [9,10,13]';
validIndices = all(not([1:length(responseCurves)]==invalidCurves),1);
responseCurves = responseCurves(validIndices);

% Execute ARCGen. Do not output processed data structure. 
[charAvg, innCorr, outCorr, ~] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts);

figure('Name','Interspinous Ligament');
title('Interspinous Ligament')
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Construct Mattucci Average. Parameters from Mattucci & Cronin (2015)
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

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pMattucci,pAvg,pCorr], 'Location', 'Best')
xlim([0,12])
ylim([0,150])
grid on
xlabel('Deflection (mm)')
ylabel('Force (N)')