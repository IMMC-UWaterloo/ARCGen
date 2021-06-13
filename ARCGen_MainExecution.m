%% ARCGen Main Execution File
% 
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           13-Jun-2021
% Updated By:     
% Date:           
% Version:        MATLAB R2020b (older versions not guaranteed)
%
% TODO: add description and syntax

%% Initialization
fclose all;
close all;
clear;
clc;

% Load data
load('Data/Lessley Parabolas/Lessley_Parabola_Processed.mat')
invalidCurves = [4];

% Remove invalid curves
validIndices = not([1:length(responseCurves)]==invalidCurves);

% Run ARCGen using ellipsoids
charAvg = ARCGen_Ellipsoidal(responseCurves(validIndices),...
    'Diagnostics','on',...
    'nResamplePoints',100);