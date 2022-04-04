%% Kroell Thoracic Impact Test Case
% This test cases demonstrates the effacacy of ARCGen on hysteretic or 
% load-unload signals. In this case, the number of warping potns must be
% chosen by examining the (x, s^hat) and (y, s^hat) for critical features
% and inflection points. 
%
% Please refer to test case readme document for more details. 
%
% Dataset Citation:
%    Kroell, C. K., Schneider, D. C., & Nahum, A. M. (1971). "Impact 
%        Tolerance and Response of the Human Thorax." SAE Technical Papers.
%        https://doi.org/10.4271/710851</div>
%
%    Lobdell, T. E., Kroell, C. K., Schneider, D. C., Hering, W. E., & 
%        Hahum, A. M. (1972). "Impact Response of the Human Thorax." In 
%        W. F. King & H. J. Mertz (Eds.), "Human Impact Response: 
%       Measurement and Simulation" (pp. 201–245). Springer Science + 
%       Business Media.
%
% Copyright (c) 2022 Devon C. Hartlen 

%% MATLAB initialiazation
fclose all;
close all;
clear;
clc;

% Set resample points and corridor resolution. 500 points will demonstrate
% runtime difference between `UseParallel` set to `on` and `off` without
% taking excessively long to run. 
nResample = 500;
nCorrPts = 500;

addpath('../') % Ensure ARCGen is on execution path

%% Load Lobdell corridors. 
% Lobdell corridors are presented in the literature as skeletal deflection
% and force response, not the measured values. Corridors are back to total
% force and defleciton here. 
lobdellCorrs = readmatrix('Kroell 1971 Thorax Response/Lobdell 16mph Corridors.csv');
charAvgLobdell = lobdellCorrs(:,1:2);
charAvgLobdell = [charAvgLobdell(:,1)+0.5, charAvgLobdell(:,2)-150];
innCorrLobdell = lobdellCorrs(~isnan(lobdellCorrs(:,3)),3:4);
innCorrLobdell = [innCorrLobdell(:,1)+0.5, innCorrLobdell(:,2)-150];
outCorrLobdell = lobdellCorrs(~isnan(lobdellCorrs(:,5)),5:6);
outCorrLobdell = [outCorrLobdell(:,1)+0.5, outCorrLobdell(:,2)-150];

% Plot the Kroell signals on thier own without corridors. 
load('Kroell 1971 Thorax Response/KroellThoraxResponse_1971.mat')
figure('Name','No Normalization');
hold on;
cmap = cbrewer2('Paired',length(responseCurves));
for iPlot = 1:length(responseCurves)
    plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName',responseCurves(iPlot).specId,...
        'LineWidth',1,'Color',cmap(iPlot,:));
end

legend('Location', 'Best')
xlim([0, 4.5])
ylim([0, 1200])
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')


%% Generate Corridors
% Generate corridors. Use three warping control points with a warping
% penalty of 1e-2. Diagonstics are turned on to see how x and y evolve with
% normalized arc-length
[charAvgNorm, innCorrNorm, outCorrNorm, proCurveDataNorm] = ...
    arcgen(responseCurves,...
    'Diagnostics', 'on', ...
    'nResamplePoints', nResample,...
    'CorridorRes', nCorrPts,...
    'nWarpCtrlPts', 3,...
    'warpingPenalty', 1e-2,...
    'UseParallel', 'on');

% Plot input signals with ARCGen and Lobdell corridors. 
figure('Name','Normalization');
hold on;
for iPlot = 1:length(responseCurves)
    pExp = plot(responseCurves(iPlot).data(:,1),...
        responseCurves(iPlot).data(:,2),...
        'DisplayName','Exp.',...
        'LineWidth',1,'Color',0.7.*[1,1,1]);
end

% Plot Lobdell average and corridors
pAvgLob = plot(charAvgLobdell(:,1),charAvgLobdell(:,2),'-',...
    'DisplayName','Char. Avg. - Lobdell','MarkerSize',16,...
    'LineWidth',2.0,'Color',[55,126,184]./255);
pCorrLob = plot(innCorrLobdell(:,1),innCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Corridors - Lobdell',...
    'LineWidth',1.0,'Color',[55,126,184]./255);
p = plot(outCorrLobdell(:,1),outCorrLobdell(:,2),'o-','MarkerSize',6,...
    'DisplayName','Outer',...
    'LineWidth',1.0,'Color',[55,126,184]./255);

% Plot ARCGen average and corridors. 
pAvg = plot(charAvgNorm(:,1),charAvgNorm(:,2),'.-',...
    'DisplayName','Char. Avg.','MarkerSize',16,...
    'LineWidth',2.5,'Color',[0,0,0]);
pCorr = plot(innCorrNorm(:,1),innCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Corridors',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);
p = plot(outCorrNorm(:,1),outCorrNorm(:,2),'.-','MarkerSize',16,...
    'DisplayName','Outer',...
    'LineWidth',2.0,'Color',[255, 213, 79]./255);

legend([pExp,pAvg,pCorr,pAvgLob,pCorrLob], 'Location', 'Best')
xlim([0, 4.5])
ylim([0, 1200])
grid on
xlabel('Deflection (in)')
ylabel('Force (lb)')