%% Pre-process Input Signals
%
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           13-Jun-2021
% Updated By:     D.C. Hartlen, M.ASc, EIT     
% Date:           18-Sep-2021           
% Version:        MATLAB R2020b (older versions not guaranteed)
%
% This script pre-processes input signals which will be subsequently used
% for to generate a characteristic average and response corridors.  
% Preprocessing amalgoamtes several curves, ensures the validity of said
% curves, and saves the data into a single MATLAB data file for later use. 
%
% Input signals must be saved in individual CSV files, with data saved in
% columns. The user can set which two column indices are used as the input
% x and y data. 
%
% Corridor generation is performed in a separate script
%
% Corridor generation scripts requires that input data be organized using a
% structure array. The structure array must have two entries per response
% curve
%   + data: [n,2] array of x-y data
%   + specId: A character string used as a specimen identifier. 
%
% This script has four options to specify "specId" in a programatic
% fashion. This is defined using "flagAlterSpecID".
%   + "No": "specId" is taken directly from the file name of the .csv
%   + "RemoveUnderscore": "specId" is the file name of the .csv with
%        underscores replaced with spaces
%   + "Squential": "specId" is defined sequential with "sequentialBase"
%        used as a prefix
%   + "Manual": "specId" is defined using the cell array "manualSpecIds".
%        "manualSpecId" must be the same length as .csv file being 
%        processed. 

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
inputSignals = struct([]); % initialization
for iFile = 1:length(inputFilenames)
    curveData = readmatrix(inputFilenames(iFile).name); % R2020a required
    inputSignals(iFile).specId = inputFilenames(iFile).specId;
    inputSignals(iFile).data = curveData(:,indicesCurves);
end
    
%% Save response curves to file
uisave({'inputSignals'})

