%% NBDL Head Kinematics, 15g Frontal Decelerations, Head Displacements
% This script was used to generate plots for Hartlen & Cronin (2022). This
% script also serves as a great example of how the number of warping points
% and warping penalty factor influences the resulting corridors. The
% parameters provided with this script were found to produce the best
% results. 
%
% This script produces corridors for head displacements. Please refer to 
% "TestCase_NBDL_15gFrontal_Disp.m" for head accelerations. 
%
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

%% Plot Head X-Displacement
% Load the point-wise average and standard deviation of timeseries data
load('NBDL 15g Frontal/NBDL 15g Frontal - Head XDisp_PointWiseAvg.mat');
nbdlAvg = pointAvg;
nbdlInner = innerCorr;
nbdlOuter = outerCorr;

% Load input signals
load('NBDL 15g Frontal/NBDL 15g Frontal - Head XDisp.mat')

% Execute ARCGen with 3 warping points and a 1e-2 penalty factor. Use
% parallel computing for faster execution. 
[charAvg, innCorr, outCorr, proCurveData] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts,...
    'nWarpCtrlPts', 3,...
    'warpingPenalty', 1e-2,...
    'UseParallel', 'on');

% Plot signals, pointwise average, and ARCGen results. 

figure('Name','Head X-Displacement');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(proCurveData(iPlot).data(:,1),...
        proCurveData(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Plot pointwise average and standard deviation corridor.
pAvgNbdl = plot(nbdlAvg(:,1),nbdlAvg(:,2),'-',...
    'DisplayName','Pointwise Avg.','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrNbdl = plot(nbdlInner(:,1),nbdlInner(:,2),'-','MarkerSize',6,...
    'DisplayName','Pointwise Corridor',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(nbdlOuter(:,1),nbdlOuter(:,2),'-','MarkerSize',6,...
    'DisplayName','Pointwise Corridor',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

% Plot ARCGen average and corridors
pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

title('Head X Displacement')
grid on
xlabel('Time (ms)')
ylabel('Displacement (mm)')

%% Plot Head Z-Displacement
% Load the point-wise average and standard deviation of timeseries data
load('NBDL 15g Frontal/NBDL 15g Frontal - Head ZDisp_PointWiseAvg.mat');
nbdlAvg = pointAvg;
nbdlInner = innerCorr;
nbdlOuter = outerCorr;

% Load input signals
load('NBDL 15g Frontal/NBDL 15g Frontal - Head ZDisp.mat')

% Execute ARCGen with 2 warping points and a 1e-2 penalty factor. Use
% parallel computing for faster execution. 
[charAvg, innCorr, outCorr, proCurveData] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts,...
    'nWarpCtrlPts', 2,...
    'warpingPenalty', 1e-2,...
    'UseParallel', 'on');

% Plot signals, pointwise average, and ARCGen results. 
figure('Name','Head Z-Displacement');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(proCurveData(iPlot).data(:,1),...
        proCurveData(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Plot pointwise average and standard deviation corridor. 
pAvgNbdl = plot(nbdlAvg(:,1),nbdlAvg(:,2),'-',...
    'DisplayName','Pointwise Avg.','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrNbdl = plot(nbdlInner(:,1),nbdlInner(:,2),'-','MarkerSize',6,...
    'DisplayName','Pointwise Corridor',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(nbdlOuter(:,1),nbdlOuter(:,2),'-','MarkerSize',6,...
    'DisplayName','Pointwise Corridor',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

% Plot ARCGen average and corridors
pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

title('Head Z Displacement')
grid on
xlabel('Time (ms)')
ylabel('Displacement (mm)')

%% Plot Head Y Angular Displacement
% Load the point-wise average and standard deviation of timeseries data
load('NBDL 15g Frontal/NBDL 15g Frontal - Head YrotDisp_PointWiseAvg.mat');
nbdlAvg = pointAvg;
nbdlInner = innerCorr;
nbdlOuter = outerCorr;

% Load input signals
load('NBDL 15g Frontal/NBDL 15g Frontal - Head YrotDisp.mat')

% Execute ARCGen with 3 warping points and a 1e-2 penalty factor. Use
% parallel computing for faster execution. 
[charAvg, innCorr, outCorr, proCurveData] = ...
    arcgen(responseCurves,...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts,...
    'nWarpCtrlPts', 3,...
    'warpingPenalty', 1e-2,...
    'UseParallel', 'on');

% Plot signals, pointwise average, and ARCGen results. 
figure('Name','Head Y Angular Displacement');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(proCurveData(iPlot).data(:,1),...
        proCurveData(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Plot pointwise average and standard deviation corridor. 

pAvgNbdl = plot(nbdlAvg(:,1),nbdlAvg(:,2),'-',...
    'DisplayName','Pointwise Avg.','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrNbdl = plot(nbdlInner(:,1),nbdlInner(:,2),'-','MarkerSize',6,...
    'DisplayName','Pointwise Corridor',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(nbdlOuter(:,1),nbdlOuter(:,2),'-','MarkerSize',6,...
    'DisplayName','Pointwise Corridor',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

% Plot ARCGen average and corridors
pAvg = plot(charAvg(:,1),charAvg(:,2),'-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorr(:,1),innCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorr(:,1),outCorr(:,2),'-','MarkerSize',16,...
    'DisplayName','Corridor',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

title('Head Y Rotation')
grid on
xlabel('Time (ms)')
ylabel('Angular Displacement (deg)')

