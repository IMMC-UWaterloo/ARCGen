%% ARCGen - Arc-length Response Corridor Generator - Ellisoid-Based
% 
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           27-Jun-2021
% Updated By:     
% Date:           
% Version:        MATLAB R2020b (older versions not guaranteed)
%
% TODO: add description and syntax

function [charAvg, innerCorr, outerCorr, processedCurveData] = ...
    ARCGen_Ellipsoidal(responseCurves,varargin)

%% Setup Name-Value Argument parser
nvArgObj = inputParser;
addParameter(nvArgObj, 'nResamplePoints',   100);
addParameter(nvArgObj, 'Diagnostics',       'off');
addParameter(nvArgObj, 'InvalidCurves',     []);
addParameter(nvArgObj, 'CorridorScaleFact', 1);
addParameter(nvArgObj, 'NormalizeCurves',   'off');
addParameter(nvArgObj, 'HandleOutliers',    'off');
addParameter(nvArgObj, 'DeviationFact',     2);
addParameter(nvArgObj, 'EllipseKFact',      1);
addParameter(nvArgObj, 'CorridorRes',       100);
nvArgObj.KeepUnmatched = true;
parse(nvArgObj,varargin{:});

nvArg = nvArgObj.Results;  % Structure created for convenience

%% Compute arclength based on input curve datapoints
% Do not perform normalization
if strcmp(nvArg.NormalizeCurves,'off')
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

% Perform curve normalization
else
    % Extract max of x and y data
    for iCurve = 1:length(responseCurves)
        tempMax = max(responseCurves(iCurve).data,[],1);
        responseCurves(iCurve).xMax = tempMax(1);
        responseCurves(iCurve).yMax = tempMax(2);
    end
    % Decision: group mean? group max? I think mean. 
    xNorm = mean([responseCurves.xMax]);
    yNorm = mean([responseCurves.yMax]);
    % Normalize the axis of each curve, then do arc-length calcs
    for iCurve = 1:length(responseCurves)
        temp = responseCurves(iCurve).data; % Temporary for conveinence
        % normalize by simple division
        temp = [temp(:,1)./xNorm, temp(:,2)./yNorm];
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
        responseCurves(iCurve).xNormMax = tempMax(1);
        responseCurves(iCurve).yNormMax = tempMax(2);
        % Remove spurious duplicates
        [~,index,~] = unique(responseCurves(iCurve).data(:,4));
        responseCurves(iCurve).data = responseCurves(iCurve).data(index,:);
    end      
end

% Compute mean and median arc-length deviation
meanAlen = mean([responseCurves.maxAlen]);
for iCurve=1:length(responseCurves)
    responseCurves(iCurve).meanDevs = ...
        responseCurves(iCurve).maxAlen-meanAlen;
end

medianAlen = median([responseCurves.maxAlen]);
for iCurve=1:length(responseCurves)
    responseCurves(iCurve).medianDev = ...
        responseCurves(iCurve).maxAlen-medianAlen;
end

%% Begin handling of outliers
switch nvArg.HandleOutliers
    case 'RemoveExtraneous'
        % Use median absolute deviation
        indexInvalid = abs([responseCurves.medianDev]) > ...
            nvArg.DeviationFact*median(abs([responseCurves.medianDev]));
        responseCurves(indexInvalid) = [];
        
    case 'CropToShortest'
        minAlen = min([responseCurves.maxAlen]);
        % Crop each curve to shortest and fix alen normalization
        for iCurve = 1:length(responseCurves)
            index = responseCurves(iCurve).data(:,3) <= minAlen;
            responseCurves(iCurve).data = ...
                responseCurves(iCurve).data(index,:);
            
            temp = responseCurves(iCurve).data; % Temporary for conveinence
            % normalize by simple division
            temp = [temp(:,1)./xNorm, temp(:,2)./yNorm];
            % Compute arc-length between each data point
            segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
                + (temp(1:end-1,2)-temp(2:end,2)).^2);
            alen = cumsum([0;segments]);
            
            responseCurves(iCurve).data(:,3) = alen;
            % Compute normalized arc-length
            responseCurves(iCurve).maxAlen = max(alen);
            responseCurves(iCurve).data(:,4) = ...
                alen./responseCurves(iCurve).maxAlen;
            % Determine max [x,y] data
            tempMax = max(temp,[],1);
            responseCurves(iCurve).xNormMax = tempMax(1);
            responseCurves(iCurve).yNormMax = tempMax(2);
            % Remove spurious duplicates
            [~,index,~] = unique(responseCurves(iCurve).data(:,4));
            responseCurves(iCurve).data = responseCurves(iCurve).data(index,:);
        end
    
    case 'CropToDeviationFactor'
        lenCropAlen = median([responseCurves.maxAlen]) + ...
            median(abs([responseCurves.medianDev])).*nvArg.DeviationFact;
        for iCurve = 1:length(responseCurves)
            index = responseCurves(iCurve).data(:,3) <= lenCropAlen;
            responseCurves(iCurve).data = ...
                responseCurves(iCurve).data(index,:);
            
            temp = responseCurves(iCurve).data; % Temporary for conveinence
            % normalize by simple division
            temp = [temp(:,1)./xNorm, temp(:,2)./yNorm];
            % Compute arc-length between each data point
            segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
                + (temp(1:end-1,2)-temp(2:end,2)).^2);
            alen = cumsum([0;segments]);
            
            responseCurves(iCurve).data(:,3) = alen;
            % Compute normalized arc-length
            responseCurves(iCurve).maxAlen = max(alen);
            responseCurves(iCurve).data(:,4) = ...
                alen./responseCurves(iCurve).maxAlen;
            % Determine max [x,y] data
            tempMax = max(temp,[],1);
            responseCurves(iCurve).xNormMax = tempMax(1);
            responseCurves(iCurve).yNormMax = tempMax(2);
            % Remove spurious duplicates
            [~,index,~] = unique(responseCurves(iCurve).data(:,4));
            responseCurves(iCurve).data = responseCurves(iCurve).data(index,:);
        end
        
    % Weighted averages and standard deviations later on. Nothing here. 
    otherwise

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

switch nvArg.HandleOutliers
    case 'WeightedAverage'
        weightFact = ...
            median(abs([responseCurves.medianDev])) * nvArg.DeviationFact;
        for iCurve = 1:length(responseCurves)
            responseCurves(iCurve).weight = max(0,...
                (weightFact - abs(responseCurves(iCurve).medianDev))/weightFact);
        end
        
        for iPoints = 1:nvArg.nResamplePoints
            clear temp;
            for iCurve = 1:length(responseCurves)
                temp(iCurve,:) = responseCurves(iCurve).normalizedCurve(iPoints,2:3);
            end
            charAvg(iPoints,:) = sum([responseCurves.weight]'.*temp())./...
                sum([responseCurves.weight]);
            nNonZero = length(responseCurves([responseCurves.weight]'>0));
            stdevData(iPoints,:) = sqrt(...
                sum([responseCurves.weight]'.*(charAvg(iPoints,:)-temp).^2)./...
                ((nNonZero-1)/nNonZero)./sum([responseCurves.weight]));
        end
        
    % Any other outlier handling
    otherwise
        for iPoints=1:nvArg.nResamplePoints
            clear temp; % probably cleaner way to do this.
            % collect specific point from each data curve
            for iCurve=1:length(responseCurves)
                temp(iCurve,:) = responseCurves(iCurve).normalizedCurve(iPoints,2:3);
            end
            charAvg(iPoints,:) = mean(temp,1);
            stdevData(iPoints,:) = std(temp,1);
        end
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
            'color',cmap(iCurve,:),...
            'DisplayName',responseCurves(iCurve).specId)
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

%% Begin marching squares algorithm

% Create grids based on upper and lower of characteristic average plus 120%
% of maximum standard deviation
[xx,yy] = meshgrid(...
    linspace(min(charAvg(:,1)) - 2.0*max(stdevData(:,1)), ...
        max(charAvg(:,1)) + 2.0*max(stdevData(:,1)), nvArg.CorridorRes),...
    linspace(min(charAvg(:,2)) - 2.0*max(stdevData(:,2)), ...
        max(charAvg(:,2)) + 2.0*max(stdevData(:,2)), nvArg.CorridorRes));
zz = zeros(size(xx));   % initalize grid of ellipse values

% For each grid point, find the max of each standard deviation ellipse
for iPt = 1:nvArg.CorridorRes
    for jPt = 1:nvArg.CorridorRes
        for kEllip = 1:size(charAvg,1)
            zz(iPt,jPt) = max(zz(iPt,jPt),...
                ((xx(iPt,jPt) - charAvg(kEllip,1)).^2 ./ ...
                (stdevData(kEllip,1)*nvArg.EllipseKFact).^2 ...
                + (yy(iPt,jPt) - charAvg(kEllip,2)).^2 ./ ...
                (stdevData(kEllip,2)*nvArg.EllipseKFact).^2)^-1);
        end
    end
end

% The following segments is the marching squares algorith. The goal of this
% algorithm is to find the zz=1 isoline, as this represents the outer
% boundary of all elllipses. 
%
% Described in brief, this algorithm goes through each point, looking at
% its and its neighbours values. There are only 16 configurations of these
% squares or cells. Based on the configuration, add the appropriate line
% segments. This method uses linear interpolation to increase accuracy. 
lineSegments = []; % Initalize line segments

for iPt = 1:(nvArg.CorridorRes-1)  % Rows (y-axis)
    for jPt = 1:(nvArg.CorridorRes-1)   % Columns (x-axis)
        % Cell value definition
        %  1 -- 2 
        %  |    |
        %  |    |
        %  8 -- 4
        %
        % REMEMBER!!!! 
        % array(i,j) = array(rows, columns,) = array(y,x)
        
        % By carefully defining cell values and definitions, we can use
        % binary to simplify logic though a integer based switch case        % 
        cellValue = ...
            1*(zz(iPt,jPt)>1) + ...
            2*(zz(iPt+1,jPt)>1) + ...
            4*(zz(iPt+1,jPt+1)>1) + ...
            8*(zz(iPt,jPt+1)>1) + 1;
        
        switch cellValue
            case 1
                % No Vertices
            case 2
                % South-West
                lineSegments = [lineSegments;
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]];        
            case 3
                % West-North
                lineSegments = [lineSegments;
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
            case 4
                % North-South
                lineSegments = [lineSegments;
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt) ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
            case 5
                % North-East
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];        
            case 6  % Ambiguous
                % South-West
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]];
                
                % North-East
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];   
            case 7
                % West-East
                lineSegments = [lineSegments;...
                    [xx(iPt,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
            case 8
                % South - East
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];  
            case 9
                % South - East
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
            case 10
                % West-East
                lineSegments = [lineSegments;...
                    [xx(iPt,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
            case 11 % Ambiguous
                % West-North
                lineSegments = [lineSegments;
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                
                % South-East
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
            case 12
                % North-East
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];
            case 13
                % North-South
                lineSegments = [lineSegments;
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt) ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
            case 14
                % West-North
                lineSegments = [lineSegments;
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
            case 15
                % South-West
                lineSegments = [lineSegments;...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]]; 
            case 16
                % No vertices 
        end
    end
end

% After the marching squares algorithm, line segments are not sorted.
% Segments need to be sorted in order to create a proper polygon. 
% Start sorting algorithm in the "middle" of the polygon, under the
% assumption that there is fewer orphans in the middle. 
lastIndex = round(0.5*size(lineSegments,1));

% Use the index of lowest y-value to seed the envelope. Envelope is
% specifically vertices, not line segments. All vertices are duplicated in
% line segments, assuming no orphans. 
envelope = lineSegments(lastIndex,1:2);
% Add second vertex, as it uses the same line segment
envelope = [envelope; lineSegments(lastIndex,3:4)];
indexUsed = zeros(size(lineSegments,1),1);
indexUsed(lastIndex) = 1;

exitFlag = 0;   % Set exit flag fudge
% Go though all vertices looking for the next connecting face
for iVerts = 2:size(lineSegments,1)
    % For an enclosed polygon, all lines share vertices
    % Find the all repeated vertices
    foundVert12 = find(all(ismembertol(lineSegments(:,1:2), envelope(end,:)), 2));
    foundVert34 = find(all(ismembertol(lineSegments(:,3:4), envelope(end,:)), 2));
    
    % there will only ever be two points, distributed between foundVert12
    % and foundVert34. Select the vertex which is NOT the same as the last
    % vertex. This indices a new line segments. 
    if ~isempty(foundVert12)
        for iInd = 1:length(foundVert12)
            if foundVert12(iInd) ~= lastIndex
                envelope = [envelope; lineSegments(foundVert12(iInd),3:4)];
                lastIndex = foundVert12(iInd);
                indexUsed(lastIndex) = 1;
                exitFlag = 1;
                break;
            end
        end
    end
    
    % Fudge to ensure that there is an infinate loop after the first
    % condition statement. 
    if exitFlag == 1
        exitFlag = 0;
        continue;
    end
        
    if ~isempty(foundVert34)
        for iInd = 1:length(foundVert34)
            if foundVert34(iInd) ~= lastIndex
                envelope = [envelope; lineSegments(foundVert34(iInd),1:2)];
                lastIndex = foundVert34(iInd);
                indexUsed(lastIndex) = 1;
                break;
            end
        end
    end
end

% At this point, 'envelope' has all vertex points. Ordering, clockwise or
% clock-wise, is not known, but can be accounted for. 

% To divide upper and lower corridors, we first extend the characteristic
% average. By choice, this extension is a linear extrapolation based on the
% first or last two points. But we also need to include a section of the
% char. avg. curve as the corridor may not extend to start of char. arv. 

indexLength = round(0.2*length(charAvg));
lineStart = [charAvg(1,1)-1.0*max(stdevData(:,1)), ...
    interp1(charAvg(1:2,1),charAvg(1:2,2),charAvg(1,1)-1.0*max(stdevData(:,1)),'linear','extrap');
    charAvg(1:indexLength,:)];
lineEnd =  [charAvg(end-indexLength:end,:);...
    charAvg(end,1)+1.0*max(stdevData(:,1)), ...
    interp1(charAvg(end-1:end,1),charAvg(end-1:end,2),charAvg(end,1)+1.0*max(stdevData(:,1)),'linear','extrap')];

%Find intercepts to divide line using Poly
[~,~,iIntStart] = polyxpoly(envelope(:,1),envelope(:,2),...
    lineStart(:,1),lineStart(:,2))
iIntStart = iIntStart(1);

[~,~,iIntEnd] = polyxpoly(envelope(:,1),envelope(:,2),...
    lineEnd(:,1),lineEnd(:,2))
iIntEnd = iIntEnd(1);

% To divide inner or outer corridors, first determine if polygon is clockwise
% or counter-clockwise. Then, based on which index is large, separate out
% inner and outer corridor based on which intercept index is larger. 
if ispolycw(envelope(:,1),envelope(:,2))
    if iIntStart > iIntEnd
        outerCorr = [envelope(iIntStart:end,:);envelope(1:iIntEnd,:)];
        innerCorr = envelope(iIntEnd:iIntStart,:);
    else
        outerCorr = envelope(iIntStart:iIntEnd,:);
        innerCorr = [envelope(iIntEnd:end,:);envelope(1:iIntStart,:)];
    end
else
    if iIntStart > iIntEnd
        innerCorr = [envelope(iIntStart:end,:);envelope(1:iIntEnd,:)];
        outerCorr = envelope(iIntEnd:iIntStart,:);
    else
        innerCorr = envelope(iIntStart:iIntEnd,:);
        outerCorr = [envelope(iIntEnd:end,:);envelope(1:iIntStart,:)];
    end
end


%% Draw Ellipses for debug
if strcmp(nvArg.Diagnostics,'on')
    figure(); hold on;
    % Scatter plot for debug
    cmap = cbrewer2('set2',2);
    colormap(cmap);
    scatter(xx(:),yy(:),20,zz(:)>=1,'Filled');

    % plot ellipses based on standard deviation
    for iPoint=1:nvArg.nResamplePoints
        ellipse(stdevData(iPoint,1).*nvArg.EllipseKFact,...
            stdevData(iPoint,2).*nvArg.EllipseKFact,0,...
            charAvg(iPoint,1), charAvg(iPoint,2),...
            0.8.*[1,1,1]);
    end
    plot(envelope(:,1),envelope(:,2),'x-b','LineWidth',2.0)
    cmap = lines(length(responseCurves));
    for iCurve=1:length(responseCurves)
        plot(responseCurves(iCurve).data(:,1),...
            responseCurves(iCurve).data(:,2),'-',...
            'DisplayName',responseCurves(iCurve).specId,...
            'Color', cmap(iCurve,:))
    end
%     plot(charAvg(:,1),charAvg(:,2),'.-k','DisplayName','Char Avg',...
%         'LineWidth',2.0,'MarkerSize',16)
    plot(lineStart(:,1),lineStart(:,2),'.-k','DisplayName','Char Avg',...
        'LineWidth',2.0,'MarkerSize',16)
    plot(lineEnd(:,1),lineEnd(:,2),'.-k','DisplayName','Char Avg',...
        'LineWidth',2.0,'MarkerSize',16)
    
end

processedCurveData = responseCurves;
end  % End main function

% helper function to perform linear interpolation to an isovalue of 1 only
function val = interpVal(x1, y1, x2, y2)
    val = x1+(x2-x1)*(1-y1)/(y2-y1);
end
   