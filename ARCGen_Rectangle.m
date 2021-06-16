%% ARCGen - Arc-length Response Corridor Generator - Ellisoid-Based
% 
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           13-Jun-2021
% Updated By:     
% Date:           
% Version:        MATLAB R2020b (older versions not guaranteed)
%
% TODO: add description and syntax

function [charAvg, innerCorr, outerCorr, quadrants] = ARCGen_Rectangle(responseCurves,varargin)

%% Setup Name-Value Argument parser
nvArgObj = inputParser;
addParameter(nvArgObj, 'nResamplePoints',   100);
addParameter(nvArgObj, 'Diagnostics',       'off');
addParameter(nvArgObj, 'invalidCurves',     []);
addParameter(nvArgObj, 'CorridorScaleFact', 1);
addParameter(nvArgObj, 'Normalization',     'off');
nvArgObj.KeepUnmatched = true;
parse(nvArgObj,varargin{:});

nvArg = nvArgObj.Results;  % Structure created for convenience


%% Compute arclength based on input curve datapoints
for iCurve = 1:length(responseCurves)
    temp = responseCurves(iCurve).data; % Temporary for conveinence
    % Compute arc-length between each data point
    segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
        + (temp(1:end-1,2)-temp(2:end,2)).^2);
    alen = cumsum([0;segments]);
    % Append cumulative arc length to data array
    responseCurves(iCurve).data = [responseCurves(iCurve).data,alen];
    % Compute normalized arc-length
    responseCurves(iCurve).maxAlen = max(alen);
    responseCurves(iCurve).data = [responseCurves(iCurve).data,...
        alen./responseCurves(iCurve).maxAlen];
    % Determine max [x,y] data
    tempMax = max(temp,[],1);
    responseCurves(iCurve).xMax = tempMax(1);
    responseCurves(iCurve).yMax = tempMax(2);
    % Remove spurious duplicates
     [~,index,~] = unique(responseCurves(iCurve).data(:,4));
     responseCurves(iCurve).data = responseCurves(iCurve).data(index,:);
end

%% Resample response curve based on normalized arc-length
for iCurve=1:length(responseCurves)
    % Linear-interpolation for x,y data against arc-length
	cfitx = fit(responseCurves(iCurve).data(:,4),...
        responseCurves(iCurve).data(:,1),'linearinterp');
	cfity = fit(responseCurves(iCurve).data(:,4),...
        responseCurves(iCurve).data(:,2),'linearinterp');
    normAlen = linspace(0,1,nvArg.nResamplePoints)';
    % Resulting array is normalized arc-length, resampled x, resam. y
	responseCurves(iCurve).normalizedCurve = ...
        [normAlen,cfitx(normAlen),cfity(normAlen)];
end
    
%% For each resampled point, determine average and standard deviation across curves
% Initialize arrays
charAvg = zeros(nvArg.nResamplePoints,2);
stdevData = zeros(nvArg.nResamplePoints,2);

for iPoints=1:nvArg.nResamplePoints
    clear temp; % probably cleaner way to do this. 
    % collect specific point from each data curve
    for iCurve=1:length(responseCurves)
        temp(iCurve,:) = responseCurves(iCurve).normalizedCurve(iPoints,2:3);
    end
    charAvg(iPoints,:) = mean(temp,1);
    stdevData(iPoints,:) = std(temp,1);
end

%% Diagnostic: Plot normalized curves and St. Devs. 
if strcmp(nvArg.Diagnostics,'on')
    figure('Name','Diagnostic Curves');
    cmap = lines(length(responseCurves));
    % Plot normalized x,y data
    subplot(2,2,[1,2]); hold on;
    for iCurve=1:length(responseCurves)
        plot(responseCurves(iCurve).normalizedCurve(:,2),...
            responseCurves(iCurve).normalizedCurve(:,3),'.-',...
            'color',cmap(iCurve,:))
    end
    xlabel('x-data')
    ylabel('y-data')
    title('Arc-length Discretized Normalized Curves')
    % Plot normalized x data against arc-length with st. dev.
    subplot(2,2,3); hold on;
    errorbar(responseCurves(1).normalizedCurve(:,1),charAvg(:,1),...
        stdevData(:,1),'color',0.5.*[1,1,1])
    cmap = lines;
    for iCurve=1:length(responseCurves)
        plot(responseCurves(iCurve).normalizedCurve(:,1),...
            responseCurves(iCurve).normalizedCurve(:,2),'.-',...
            'color',cmap(iCurve,:))
    end
    xlabel('Normalized Arc-length')
    ylabel('x-data')
    title('Average and St.Dev. of X-Data')
    % Plot normalized y data against arc-length with st. dev.
    subplot(2,2,4); hold on;
    errorbar(responseCurves(1).normalizedCurve(:,1),charAvg(:,2),...
        stdevData(:,2),'color',0.5.*[1,1,1])
    cmap = lines;
    for iCurve=1:length(responseCurves)
        plot(responseCurves(iCurve).normalizedCurve(:,1),...
            responseCurves(iCurve).normalizedCurve(:,3),'.-',...
            'color',cmap(iCurve,:))
    end
    xlabel('Normalized Arc-length')
    ylabel('y-data')
    title('Average and St.Dev. of Y-Data')
end

%% Create Corridors
% Initialize corridors
outerCorr = zeros(nvArg.nResamplePoints,2);
innerCorr = zeros(nvArg.nResamplePoints,2);
% Initialize changeover indices
quadrants = zeros(nvArg.nResamplePoints,1);
% Append mirroed point to average for slope calculations
genAvg = [-charAvg(2,:);charAvg];
% Cycle through each resampled point
for iPoint = 1:nvArg.nResamplePoints
    % Determine directions
    curveDir = genAvg(iPoint+1,:)-genAvg(iPoint,:);
    % Start building corridors
    % If dx > 0 & dy > 0, outer = UL point, Inner = LR
    if (curveDir(1) >= 0) && (curveDir(2) >= 0)
        outerCorr(iPoint,:) = [...
            charAvg(iPoint,1)-stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)+stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        innerCorr(iPoint,:) = [...
            charAvg(iPoint,1)+stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)-stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        quadrants(iPoint) = 1;
    % If dx > 0 & dy < 0, outer = UR point, Inner = LL    
    elseif (curveDir(1) >= 0) && (curveDir(2) < 0)
        outerCorr(iPoint,:) = [...
            charAvg(iPoint,1)+stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)+stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        innerCorr(iPoint,:) = [...
            charAvg(iPoint,1)-stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)-stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        quadrants(iPoint) = 2;
    % If dx < 0 & dy < 0, outer = LR point, Inner = UL
    elseif (curveDir(1) < 0) && (curveDir(2) < 0)
        outerCorr(iPoint,:) = [...
            charAvg(iPoint,1)+stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)-stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        innerCorr(iPoint,:) = [...
            charAvg(iPoint,1)-stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)+stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        quadrants(iPoint) = 3;
    % If dx < 0 & dy > 0, outer = LL point, Inner = UR
    elseif (curveDir(1) < 0) && (curveDir(2) >= 0)
        outerCorr(iPoint,:) = [...
            charAvg(iPoint,1)-stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)-stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        innerCorr(iPoint,:) = [...
            charAvg(iPoint,1)+stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)+stdevData(iPoint,2).*nvArg.CorridorScaleFact];
        quadrants = 4;
    end
end

%% Draw rectangles (Temporary, still figuring out algorithm to extract corridors)
if strcmp(nvArg.Diagnostics,'on')
    figure(); hold on;
    % plot rectangles based on standard deviation
    for iPoint=1:nvArg.nResamplePoints
        rectangle('Position',...
            [charAvg(iPoint,1)-stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            charAvg(iPoint,2)-stdevData(iPoint,2).*nvArg.CorridorScaleFact,...
            2.*stdevData(iPoint,1).*nvArg.CorridorScaleFact,...
            2.*stdevData(iPoint,2).*nvArg.CorridorScaleFact],...
            'EdgeColor',0.85*[1,1,1])
    end
    cmap = lines;
    for iCurve=1:length(responseCurves)
        plot(responseCurves(iCurve).data(:,1),...
            responseCurves(iCurve).data(:,2),'-',...
            'DisplayName',responseCurves(iCurve).specId,...
            'Color', cmap(iCurve,:))
    end
    plot(charAvg(:,1),charAvg(:,2),'.-k','DisplayName','Char Avg',...
        'LineWidth',2.0,'MarkerSize',16)
    plot(outerCorr(:,1),outerCorr(:,2),'.-','DisplayName','Outer Corr',...
        'LineWidth',2.0,'MarkerSize',16,'Color',[89,60,31]./255)
    plot(innerCorr(:,1),innerCorr(:,2),'.-','DisplayName','Inner Corr',...
        'LineWidth',2.0,'MarkerSize',16,'Color',[89,60,31]./255)
end
   