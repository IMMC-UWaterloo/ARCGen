%% Lessley Parabolas Test Case
% This test case serves as a basic proof of concept for ARCGen, as well as
% a tutorial for new users. This test case runs ARCGen three times with
% different options to demonstrate how functionality is effected. 
%
% Please refer to test case readme document for more details. 
%
% Dataset Citation:
%    Lessley, D., Crandall, J., Shaw, G., Kent, R., & Funk, J. (2004). "A
%       normalization technique for developing corridors from individual 
%       subject responses." SAE Technical Papers.
%       https://doi.org/10.4271/2004-01-0288
%
% Copyright (c) 2022 Devon C. Hartlen 

%% MATLAB initialization
fclose all;
close all;
clear;
clc;

addpath('../'); % ensure arcgen.m is on execution path. 

%% Call 1: Basic ARCGen input options
% Load data. Note the required input format for ARCGen
load('Lessley Parabolas/Lessley_Parabola_Processed.mat')
% remove the outlier for now (located in 4th position)
invalidCurves = [4];
validIndices = not([1:length(responseCurves)]==invalidCurves);
responseCurves = responseCurves(validIndices);

% Run ARCGen, saving the average, corridors, and processed data structure
[charAvg, innCorr, outCorr,proCurveData] = ...
    arcgen(responseCurves,...
    'nResamplePoints', 200,... % How many points used to resample signals
    'CorridorRes', 100); % used to define a 100 x 100 grid for corridor extraction


% Plot input signals, characteristic average and corridors. 
figure('Name','Call 1: Basic input arguments');
hold on;
for iPlot = 1:length(proCurveData)
    pExp = plot(proCurveData(iPlot).data(:,1),...
        proCurveData(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pInn = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
pOut = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pAvg,pInn,pOut], 'Location', 'Best')
xlim([0,2.0])
ylim([0,160])
grid on
xlabel('X-Comp')
ylabel('Y-Comp')

%% Call 2: Running ARCGen without signal normalization. 
% Running without normalizing signals is strongly discouraged. 

% Load data. Note the required input format for ARCGen
load('Lessley Parabolas/Lessley_Parabola_Processed.mat')
% remove the outlier for now (located in 4th position)
invalidCurves = [4];
validIndices = not([1:length(responseCurves)]==invalidCurves);
responseCurves = responseCurves(validIndices);

% Run ARCGen, saving the average, corridors, and processed data structure
[charAvg, innCorr, outCorr,proCurveData] = ...
    arcgen(responseCurves,...
    'nResamplePoints', 200,... % How many points used to resample signals
    'CorridorRes', 100,... % used to define a 100 x 100 grid for corridor extraction
    'NormalizeSignals', 'off'); % Default: 'on'. Default highly recommended

% Plot input signals, characteristic average and corridors. 
figure('Name','Call 2: Effect of disabling mangitude scaling');
hold on;
for iPlot = 1:length(proCurveData)
    pExp = plot(proCurveData(iPlot).data(:,1),...
        proCurveData(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pInn = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
pOut = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pAvg,pInn,pOut], 'Location', 'Best')
xlim([0,2.0])
ylim([0,160])
grid on
xlabel('X-Comp')
ylabel('Y-Comp')


%% Call 3: Effect of Outlier, Demonstration of Diagontics
% An additional figure is created containing information detailing how
% re-parmeterization was accomplished. 

% Load data, do not exclude outlier
load('Lessley Parabolas/Lessley_Parabola_Processed.mat')

% Run ARCGen, saving the average, corridors, and processed data structure
[charAvg, innCorr, outCorr,proCurveData] = ...
    arcgen(responseCurves,...
    'nResamplePoints', 200,... % How many points used to resample signals
    'CorridorRes', 100,... % used to define a 100 x 100 grid for corridor extraction
    'Diagnostics', 'on'); % Create a diagonstic figure. 

% Plot input signals, characteristic average and corridors. 
figure('Name','Call 3: Effect of outlier');
hold on;
for iPlot = 1:length(proCurveData)
    pExp = plot(proCurveData(iPlot).data(:,1),...
        proCurveData(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

pAvg = plot(charAvg(:,1),charAvg(:,2),'.-',...
    'DisplayName','Char. Avg. - ARCGen','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pInn = plot(innCorr(:,1),innCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors - ARCGen',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
pOut = plot(outCorr(:,1),outCorr(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[196, 147, 0]./255);

legend([pExp,pAvg,pInn,pOut], 'Location', 'Best')
xlim([0, 2.5])
ylim([0, 480])
grid on
xlabel('X-Comp')
ylabel('Y-Comp')