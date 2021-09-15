%% Pre-process Input Response Curves
% 
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           13-Jun-2021
% Updated By:     D.C. Hartlen, M.ASc, EIT     
% Date:           13-Jul-2021           
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
% Actual corridor generation is performed in a separate script
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

addpath('../../ThirdPartyFunctions') % Path to 3rd party functions

sheetNames = {...
    'X-Accel',...
    'Z-Accel',...
    'Yrot-Accel',...
    'X-Disp2',...
    'Z-Disp5',...
    'Yrot-Disp8',...
    'T1_X-Accel47',...
    'T1_Z-Accel49',...
    'T1_X-Disp30',...
    'T1_Z-Disp33',...
    'T1 - Yrot Angle80 CorrectSpline',...
    };

idNames = {...
    'Head XAccel',...
    'Head ZAccel',...
    'Head YrotAccel',...
    'Head XDisp',...
    'Head ZDisp',...
    'Head YrotDisp',...
    'T1 XAccel',...
    'T1 ZAccel',...
    'T1 XDisp',...
    'T1 ZDisp',...
    'T1 YrotDisp',...
    };

for iSheet = 1:length(sheetNames)

data = readmatrix('8g Frontal.xlsx','sheet',sheetNames{iSheet});

specIdCols = [7:16];
timeCol = 2;
avgCol = 3;
stDevCol = 4;

%% Load response curves
figure(); hold on;
% Specify x,y columns of datafile to be loaded
responseCurves = struct([]); % initialization
for iCurve = 1:length(specIdCols)
    responseCurves(iCurve).specId = num2str(data(1,specIdCols(iCurve)));
    responseCurves(iCurve).data = [data(2:end,timeCol),...
        data(2:end,specIdCols(iCurve))];
    responseCurves(iCurve).data(isnan(responseCurves(iCurve).data)) = 0;
    plot(responseCurves(iCurve).data(:,1),responseCurves(iCurve).data(:,2),...
        '.-','DisplayName',responseCurves(iCurve).specId)
end

pointAvg = [data(2:end,timeCol),data(2:end,avgCol)];
pointAvg(isnan(pointAvg)) = 0;
innerCorr = [data(2:end,timeCol),data(2:end,avgCol)-data(2:end,stDevCol)];
innerCorr(isnan(innerCorr)) = 0;
outerCorr = [data(2:end,timeCol),data(2:end,avgCol)+data(2:end,stDevCol)];
outerCorr(isnan(outerCorr)) = 0;

plot(pointAvg(:,1),pointAvg(:,2),'DisplayName','Pointwise Avg',...
    'Color',[0,0,0],'LineWidth',2.0)
plot(innerCorr(:,1),innerCorr(:,2),'DisplayName','Inner Corridor',...
    'Color',0.8.*[1,1,1],'LineWidth',2.0)
plot(outerCorr(:,1),outerCorr(:,2),'DisplayName','Outer Corridor',...
    'Color',0.8.*[1,1,1],'LineWidth',2.0)


%% Filter curves
% TODO: Add filtering

%% Save response curves to file
save(['NBDL 8g Frontal - ' idNames{iSheet} '.mat'], 'responseCurves')
save(['NBDL 8g Frontal - ' idNames{iSheet} '_PointWiseAvg.mat'],...
    'pointAvg','innerCorr','outerCorr')
end