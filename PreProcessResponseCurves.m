%% Pre-process Input Response Curves
% 
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           13-Jun-2021
% Updated By:     
% Date:           
% Version:        MATLAB R2020b (older versions not guaranteed)
%
% This script pre-processes input curves which will be subsequently used
% for to generate a characteristic average and response corridors. This
% preprocessing amalgoamtes several curves, ensures the validity of said
% curves, optionally filters the curves, and saves the data into a single
% MATLAB data file for later use. 
%
% Input curves must be saved in individual CSV files, with data saved in
% columns. The user can set which two column indices are used as the input
% x and y data. 
%
% Acutal response corridor fitting is accomplished in a separate script

%% Initialization
fclose all;
close all;
clear;
clc;

addpath('ThirdPartyFunctions') % Path to 3rd party functions

%% Select desired data files to be processed
inputFilenames = uipickfiles('Output','struct');
% correct file names
for iFile = 1:length(inputFilenames)
    [~,name,ext] = fileparts(inputFilenames(iFile).name);
    inputFilenames(iFile).specId = name;
end

%% Alter specimen ID from file name (if desired)
% Cases: 'No','RemoveUnderscore','Sequential','Manual'
flagAlterSpecID = 'RemoveUnderscore'; 
sequentialBase = 'ID';
manualSpecIds = {...
    'ID1';...
    'ID2';...
    };
switch flagAlterSpecID
    case 'No'
        disp('Skipping rename')
    case 'RemoveUnderscore'
        for iFile = 1:length(inputFilenames)
            inputFilenames(iFile).specId = ...
                replace(inputFilenames(iFile).specId,'_(x)','');
            inputFilenames(iFile).specId = ...
                replace(inputFilenames(iFile).specId,'_',' ');
        end
    case 'Sequential'
        for iFile = 1:length(inputFilenames)
            inputFilenames(iFile).specId = ...
                [sequentialBase ' ' num2str(iFile,'%3d')];
        end
    case 'Manual'
        if length(manualSpecIds) ~= length(inputFilenames)
            error('Not enough manual IDs specified')
        else
            for iFile = 1:length(inputFilenames)
                inputFilenames(iFile).specId = manualSpecIds{iFile};
            end
        end
end

%% Load response curves
% Specify x,y columns of datafile to be loaded
indicesCurves = [1,2];
responseCurves = struct([]); % initialization
for iFile = 1:length(inputFilenames)
    curveData = readmatrix(inputFilenames(iFile).name); % R2020a required
    responseCurves(iFile).specId = inputFilenames(iFile).specId;
    responseCurves(iFile).data = curveData(:,indicesCurves);
end
    
%% Filter curves
% TODO: Add filtering

%% Save response curves to file
uisave({'responseCurves'})

