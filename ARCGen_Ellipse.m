%% ARCGen_Ellipse - Arc-length Response Corridor Generator - Ellise-Based
% 
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           27-Jun-2021
% Updated By:     D.C. Hartlen, M.ASc, EIT
% Date:           13-Jul-2021
% Version:        MATLAB R2020b (older versions not guaranteed)
%
% ARCGen, short for Arc-length Response Corridor Generation, provides
% automated calculation of a characteristic average and response corridors
% on input curves regardless of if said curves are non-monotonic or
% hystertic. This is accomplished by re-parameterizing input curves based
% on arc-length. Corridors are extracted using a marching squares
% algorithm.
%
% This function has one mandatory input, four outputs, and many optional
% inputs. Optional inputs are defined using name-value pair arguments. 
%
% Usage notes: It is common to see errors when running this function if the
% number of resampling points or corridor extraction grid is too sparse.
% This is error also occurs if standard deviation in a particular direction
% is too small for subsequent ellipses to overlap significantly. Typically, 
% this problem will result in an error pointing to a line 650+. A quick fix
% is to keep bumping up 'nResamplePoints' and 'CorridorRes' until the error
% goes array. Turning 'Diagnostics' 'on' can help identify these issues. 
%
% MANDATORY INPUTS:
% -----------------
% responseCurves: A [nCurve,2] structured array consisting of the following
%       entries. Entries are case-senstive
%   + data: an [n,2] array containing ordered x-y data
%   + specId: character array containing an identifier for each curve
% 
% OPTIONAL INPUTS:
% ----------------
% nResamplePoints: integer defining the number of points used to
%       re-parameterize input curves. Default: 100. 
% CorridorRes: integeer defining the number of grid points used for the
%       marching squares algorithm. The sampling grid for the marching
%       squares algorithm extends 120% of extreme corridors. This parameter
%       defines the number of points along each side of the grid. 
%       Default: 100. It is common to increase this significantly. 
% NormalizeCurves: character arry used to turn on curve normalization.
%       Options: 'on' (default), 'off'
% EllipseKFact: float used to scale the major and minor axis of the
%       ellipses used for corridor generation. This value corrisponds to 
%       the square root of the chi-squared CDF. Default: 1.0 (creates 
%       corridors one standard deviation along the x and y axes)
% Diagnostics: character array used to activate diagnostic plots. Useful
%       for debugging errors. Options: 'on', 'off' (default)
% MinCorridorWidth: Factor used to enforce a minimum corridor width. Any
%       st.dev. less than 'MinCorridorFactor'*max(st.dev.) is replaced with
%       'MinCorridorFactor'*max(st.dev.). x & y axes are handled
%       separately. A value of 0 disables forcing minimum width.
% HandleOutliers: character array used to impliment a variety of methods to
%       handle curves of extraneous length. This is a fairly experimental
%       feature and functionality is not fully guaranteed. Generally,
%       outliers are defined based on arc-length.
%       Options are:
%   + 'off': (default) no outlier handling is performed
%   + 'CropToShortest': All input curves are cropped to the shortest
%           arc-length. Not useful if one curve is much shorter than the
%           others. 
%   + 'RemoveExtraneous': Removes curves which are significantly different
%           in length based on median deviation. Specifiying this option 
%           also requires one to define 'DeviationFact'. Curves that are 
%           have an absolute deviation large than  'DevationFact' * 
%           (mediand deviation) are not used for calculation of average or 
%           corridors. 
%   + 'CropToDeviationFactor' Crops curves which are significantly longer
%           than the median length. Specifying this option requires one to
%           define 'DeviationFact'. Curves longer than 'DevationFact' * 
%           (mediand deviation) are croped to this value. Curves which are
%           shorter than median are not effected at all. 
%   + 'WeightedAverage': Characteristic average and corridors are computed
%           using a weighted mean and standard devation of input curves
%           based on how far the arc-length of each curve is from the
%           median values. Specifying this option also requires one to
%           define 'DeviationFact'. Curves beyond 'DeviationFact'*(median
%           arc-length) are given a weight of zero. All other curves are
%           given a weight of zero to one based on how close thier average
%           is to the median arc-length. Developer's Note: I thought this
%           was a clever solution. It doesn't work as well as I thought. 
%
% OUTPUTS:
% --------
% charAvg: an [nResamplePoints,2] array containing the computed
%       characteristic average.
% innerCorr: an [n,2] array containing points defining the inner corridor
% outerCorr: an [n,2] array containing points defining the outer corridor
% processedCurveData: a structure array that outputs processed curves and
%       some basic statistics
%
% Note on outputted corridors: corridors are not uniformly sampled due to
% the limitations of the marching squares algorithm used to extract
% corridors in an automated fashion. Additionally, it is common that the
% corridors do not extend to the start of the curves, as the standard
% devaition at thes start of curves is typically too small to provide
% sufficeint ellipse overlap for the marching squares algorithm to extract
% continuous corridors. These corridors can be resampled and extended after
% extraction. 

function [charAvg, innerCorr, outerCorr, processedCurveData] = ...
    ARCGen_Ellipse(responseCurves,varargin)

%% Setup Name-Value Argument parser
nvArgObj = inputParser;
addParameter(nvArgObj, 'nResamplePoints',   100);
addParameter(nvArgObj, 'Diagnostics',       'off');
addParameter(nvArgObj, 'InvalidCurves',     []);
addParameter(nvArgObj, 'CorridorScaleFact', 1);
addParameter(nvArgObj, 'NormalizeCurves',   'on');
addParameter(nvArgObj, 'HandleOutliers',    'off');
addParameter(nvArgObj, 'DeviationFact',     2);
addParameter(nvArgObj, 'EllipseKFact',      1);
addParameter(nvArgObj, 'CorridorRes',       100);
addParameter(nvArgObj, 'MinCorridorWidth',    0); 
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
    
elseif strcmp(nvArg.NormalizeCurves,'PositiveYPeak')
    % Find absolute largest peak in Y data
    for iCurve = 1:length(responseCurves)
        [pks,indices] = findpeaks(responseCurves(iCurve).data(:,2),...
            'MinPeakProminence', 0.05.*max(responseCurves(iCurve).data(:,2)));
        [~,peakInd] = max(abs(pks));
        responseCurves(iCurve).alignPt = ...
            [responseCurves(iCurve).data(indices(peakInd),1),...
            responseCurves(iCurve).data(indices(peakInd),2)];
        responseCurves(iCurve).alignInd = indices(peakInd);
    end
    % Debug plot for showing peaks
    figure(); hold on;
    for iCurve = 1:length(responseCurves)
        plot(responseCurves(iCurve).data(:,1),responseCurves(iCurve).data(:,2))
        plot(responseCurves(iCurve).alignPt(1),responseCurves(iCurve).alignPt(2),'kx')
        
    end
    % Calculate arc-length, then normalize [0,1] to peak, [1,2] to terminus
    % Normalize magnitudes
    for iCurve = 1:length(responseCurves)
        tempMin = min(responseCurves(iCurve).data,[],1);
        responseCurves(iCurve).xMin = tempMin(1);
        responseCurves(iCurve).yMin = tempMin(2);
        tempMax = max(responseCurves(iCurve).data,[],1);
        responseCurves(iCurve).xMax = tempMax(1);
        responseCurves(iCurve).yMax = tempMax(2);
    end
    % Decision: group mean? group max? I think mean.
    %     xNorm = mean([responseCurves.xMax]);
    %     yNorm = mean([responseCurves.yMax]);
    xBound = [mean([responseCurves.xMin]), mean([responseCurves.xMax])];
    yBound = [mean([responseCurves.yMin]), mean([responseCurves.yMax])];
    % Normalize the axis of each curve, then do arc-length calcs
    for iCurve = 1:length(responseCurves)
        temp = responseCurves(iCurve).data; % Temporary for conveinence
        % % normalize by simple division
        % temp = [temp(:,1)./xNorm, temp(:,2)./yNorm];
        % Normalize by scale and shift
        temp = [temp(:,1)./(xBound(2)-xBound(1)),...
            temp(:,2)./(yBound(2)-yBound(1))];
        % Compute arc-length between each data point
        segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
            + (temp(1:end-1,2)-temp(2:end,2)).^2);
        alen = cumsum([0;segments]);
        % Append cumulative arc length to data array
        responseCurves(iCurve).data = [responseCurves(iCurve).data,alen];
        % Compute normalized arc-length
        responseCurves(iCurve).maxAlen = max(alen);
        
        tempInd = responseCurves(iCurve).alignInd;
        scaledAlen = [...
            (alen(1:tempInd)-alen(1))./(alen(tempInd)-alen(1));
            1+(alen(tempInd+1:end)-alen(tempInd+1))./(alen(end)-alen(tempInd+1))...
            ];
        responseCurves(iCurve).data = [responseCurves(iCurve).data,...
            scaledAlen];
        % Determine max [x,y] data
        tempMax = max(abs(temp),[],1);
        responseCurves(iCurve).xNormMax = tempMax(1);
        responseCurves(iCurve).yNormMax = tempMax(2);
        % Remove spurious duplicates
        [~,index,~] = unique(responseCurves(iCurve).data(:,4));
        responseCurves(iCurve).data = responseCurves(iCurve).data(index,:);
    end

elseif strcmp(nvArg.NormalizeCurves,'NegativeYPeak')
    % Find absolute largest peak in Y data
    for iCurve = 1:length(responseCurves)
        [pks,indices] = findpeaks(-responseCurves(iCurve).data(:,2),...
            'MinPeakProminence', 0.05.*max(-responseCurves(iCurve).data(:,2)));
        [~,peakInd] = max(abs(pks));
        responseCurves(iCurve).alignPt = ...
            [responseCurves(iCurve).data(indices(peakInd),1),...
            responseCurves(iCurve).data(indices(peakInd),2)];
        responseCurves(iCurve).alignInd = indices(peakInd);
    end
    % Debug plot for showing peaks
    figure(); hold on;
    for iCurve = 1:length(responseCurves)
        plot(responseCurves(iCurve).data(:,1),responseCurves(iCurve).data(:,2))
        plot(responseCurves(iCurve).alignPt(1),responseCurves(iCurve).alignPt(2),'kx')
        
    end
    % Calculate arc-length, then normalize [0,1] to peak, [1,2] to terminus
    % Normalize magnitudes
    for iCurve = 1:length(responseCurves)
        tempMin = min(responseCurves(iCurve).data,[],1);
        responseCurves(iCurve).xMin = tempMin(1);
        responseCurves(iCurve).yMin = tempMin(2);
        tempMax = max(responseCurves(iCurve).data,[],1);
        responseCurves(iCurve).xMax = tempMax(1);
        responseCurves(iCurve).yMax = tempMax(2);
    end
    % Decision: group mean? group max? I think mean.
    %     xNorm = mean([responseCurves.xMax]);
    %     yNorm = mean([responseCurves.yMax]);
    xBound = [mean([responseCurves.xMin]), mean([responseCurves.xMax])];
    yBound = [mean([responseCurves.yMin]), mean([responseCurves.yMax])];
    % Normalize the axis of each curve, then do arc-length calcs
    for iCurve = 1:length(responseCurves)
        temp = responseCurves(iCurve).data; % Temporary for conveinence
        % % normalize by simple division
        % temp = [temp(:,1)./xNorm, temp(:,2)./yNorm];
        % Normalize by scale and shift
        temp = [temp(:,1)./(xBound(2)-xBound(1)),...
            temp(:,2)./(yBound(2)-yBound(1))];
        % Compute arc-length between each data point
        segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
            + (temp(1:end-1,2)-temp(2:end,2)).^2);
        alen = cumsum([0;segments]);
        % Append cumulative arc length to data array
        responseCurves(iCurve).data = [responseCurves(iCurve).data,alen];
        % Compute normalized arc-length
        responseCurves(iCurve).maxAlen = max(alen);
        
        tempInd = responseCurves(iCurve).alignInd;
        scaledAlen = [...
            (alen(1:tempInd)-alen(1))./(alen(tempInd)-alen(1));
            1+(alen(tempInd+1:end)-alen(tempInd+1))./(alen(end)-alen(tempInd+1))...
            ];
        responseCurves(iCurve).data = [responseCurves(iCurve).data,...
            scaledAlen];
        % Determine max [x,y] data
        tempMax = max(abs(temp),[],1);
        responseCurves(iCurve).xNormMax = tempMax(1);
        responseCurves(iCurve).yNormMax = tempMax(2);
        % Remove spurious duplicates
        [~,index,~] = unique(responseCurves(iCurve).data(:,4));
        responseCurves(iCurve).data = responseCurves(iCurve).data(index,:);
    end
    
% Perform magnitude normalization based on bounding box
elseif strcmp(nvArg.NormalizeCurves,'on')
    % Determine bounding box of individual curves
    for iCurve = 1:length(responseCurves)
        tempMin = min(responseCurves(iCurve).data,[],1);
        responseCurves(iCurve).xMin = tempMin(1);
        responseCurves(iCurve).yMin = tempMin(2);
        tempMax = max(responseCurves(iCurve).data,[],1);
        responseCurves(iCurve).xMax = tempMax(1);
        responseCurves(iCurve).yMax = tempMax(2);
    end
    xBound = [mean([responseCurves.xMin]), mean([responseCurves.xMax])];
    yBound = [mean([responseCurves.yMin]), mean([responseCurves.yMax])];
    % Normalize the axis of each curve, then do arc-length calcs
    for iCurve = 1:length(responseCurves)
        temp = responseCurves(iCurve).data; % Temporary for conveinence
        % Normalize from bounding box to [-1,1]
        temp = [temp(:,1)./(xBound(2)-xBound(1)),...
            temp(:,2)./(yBound(2)-yBound(1))];
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
        tempMax = max(abs(temp),[],1);
        responseCurves(iCurve).xNormMax = tempMax(1);
        responseCurves(iCurve).yNormMax = tempMax(2);
        % Remove spurious duplicates
        [~,index,~] = unique(responseCurves(iCurve).data(:,4));
        responseCurves(iCurve).data = responseCurves(iCurve).data(index,:);
    end
else
    error('Input Normalization Method Not Recognized')
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
            tempMax = max(abs(temp),[],1);
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
            tempMax = max(abs(temp),[],1);
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
%     indexNan = ~isnan(responseCurves(iCurve).data(:,4));
%     responseCurves(iCurve).data = responseCurves(iCurve).data(indexNan,:);
    % Linear-interpolation for x,y data against arc-length
    normAlen = linspace(0,responseCurves(iCurve).data(end,4),...
        nvArg.nResamplePoints)';
    resampX = interp1(responseCurves(iCurve).data(:,4),...
        responseCurves(iCurve).data(:,1), normAlen);
    resampY = interp1(responseCurves(iCurve).data(:,4),...
        responseCurves(iCurve).data(:,2), normAlen);
    % Resulting array is normalized arc-length, resampled x, resam. y
    responseCurves(iCurve).normalizedCurve = [normAlen, resampX, resampY];
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

%% Clamp minimum corridor width. Disabled if 'MinCorridorWidth' == 0
% Include influence of corridor scaling factor 'EllipseKFact'
if nvArg.MinCorridorWidth > 0
    % Replace any stDevData below maximum st.dev. * 'MinCorridorWidth'
    index = stdevData <...
        (nvArg.MinCorridorWidth .* max(stdevData) .* nvArg.EllipseKFact);
    stdevData(index(:,1),1) = (nvArg.MinCorridorWidth .* nvArg.EllipseKFact...
        .* max(stdevData(:,1)));
    stdevData(index(:,2),2) = (nvArg.MinCorridorWidth .* nvArg.EllipseKFact...
        .* max(stdevData(:,2)));
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
    legend()
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
    
    % Plot ellipses    
    figure('Name','Overlaid Ellipses and Corridors'); hold on;
    % Scatter plot for debug
    cmap = cbrewer2('set2',2);
    colormap(cmap);
%     scatter(xx(:),yy(:),20,zz(:)>=1,'Filled');

    % plot ellipses based on standard deviation
    for iPoint=1:nvArg.nResamplePoints
        ellipse(stdevData(iPoint,1).*nvArg.EllipseKFact,...
            stdevData(iPoint,2).*nvArg.EllipseKFact,0,...
            charAvg(iPoint,1), charAvg(iPoint,2),...
            0.8.*[1,1,1]);
    end

    cmap = lines(length(responseCurves));
    for iCurve=1:length(responseCurves)
        plot(responseCurves(iCurve).data(:,1),...
            responseCurves(iCurve).data(:,2),'-',...
            'DisplayName',responseCurves(iCurve).specId,...
            'Color', cmap(iCurve,:))
    end
    
    plot(charAvg(:,1),charAvg(:,2),'.-k','DisplayName','Char Avg',...
        'LineWidth',2.0,'MarkerSize',16)
end

%% Begin marching squares algorithm

% Create grids based on upper and lower of characteristic average plus 120%
% of maximum standard deviation
scaleFact = 1.2*nvArg.EllipseKFact;
[xx,yy] = meshgrid(...
    linspace(min(charAvg(:,1)) - scaleFact*max(stdevData(:,1)), ...
        max(charAvg(:,1)) + scaleFact*max(stdevData(:,1)),...
        nvArg.CorridorRes),...
    linspace(min(charAvg(:,2)) - scaleFact*max(stdevData(:,2)), ...
        max(charAvg(:,2)) + scaleFact*max(stdevData(:,2)),...
        nvArg.CorridorRes));
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
                centerVal = mean([zz(iPt,jPt), zz(iPt+1,jPt), zz(iPt+1,jPt+1), zz(iPt, jPt+1)]);
                if centerVal >= 1
                    % West-North
                    lineSegments = [lineSegments;
                        [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                        interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                    % South - East
                    lineSegments = [lineSegments;...
                        [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                        xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
                else
                    % South-West
                    lineSegments = [lineSegments;...
                        [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                        xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]];
                    % North-East
                    lineSegments = [lineSegments;...
                        [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                        xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];
                end
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
                centerVal = mean([zz(iPt,jPt), zz(iPt+1,jPt), zz(iPt+1,jPt+1), zz(iPt, jPt+1)]);
                if centerVal >= 1
                    % South-West
                    lineSegments = [lineSegments;...
                        [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                        xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]];
                    % North-East
                    lineSegments = [lineSegments;...
                        [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                        xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];
                else
                    % West-North
                    lineSegments = [lineSegments;
                        [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                        interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                    % South-East
                    lineSegments = [lineSegments;...
                        [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                        xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
                end
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
%
% One issues with very small ellipses is orphan envelopes can occur. This
% sorting algorithm attempts to find all envelopes. Final corridor
% extraction is only done on the largest envelope. Largest is defined by
% the number of vertices. 
%
% Start sorting algorithm in the "middle" of the polygon, under the
% assumption that there is fewer orphans in the middle. But this algorithm
% should be fine if it is an orphan. 
iEnvelope = 1;
indexUsed = zeros(size(lineSegments,1),1);
lastIndex = round(0.5*size(lineSegments,1));

while ~all(indexUsed==1)
    % Use the index of lowest y-value to seed the envelope. Envelope is
    % specifically vertices, not line segments. All vertices are duplicated in
    % line segments, assuming no orphans.
    envelope = lineSegments(lastIndex,1:2);
    % Add second vertex, as it uses the same line segment
    envelope = [envelope; lineSegments(lastIndex,3:4)];
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
    
    envelope = unique(envelope,'stable','rows');
    individualEnvelopes{iEnvelope} = envelope;
    iEnvelope = iEnvelope+1;
    clear envelope;
    lastIndex = find(indexUsed==0,1);
    % if isempty(lastIndex)
    %     break
    % end
end

% Choose the biggest envelope, based on number of indices. 
for iEnv = 1:length(individualEnvelopes)
    envSize(iEnv) = size(individualEnvelopes{iEnv},1);
end
[~,ind] = max(envSize);
envelope = individualEnvelopes{ind};

% For debugging, plot all envelopes
if strcmp(nvArg.Diagnostics,'on')
    for iEnv = 1:length(individualEnvelopes)
            plot(individualEnvelopes{iEnv}(:,1),...
            individualEnvelopes{iEnv}(:,2),'.-b','LineWidth',1.0)   
    end
end

% At this point, 'envelope' has all vertex points. Ordering, clockwise or
% clock-wise, is not known, but can be accounted for. 

% To divide upper and lower corridors, we first extend the characteristic
% average. By choice, this extension is a linear extrapolation based on the
% first or last two points. But we also need to include a section of the
% char. avg. curve as the corridor may not extend to start of char. arv. 

aLenInterval = 1./nvArg.nResamplePoints;
indexLength = round(0.2*length(charAvg));

lineStart = [...
    interp1([0,aLenInterval],charAvg(1:2,1), -0.20,'linear','extrap'),...
    interp1([0,aLenInterval],charAvg(1:2,2), -0.20,'linear','extrap');...
    charAvg(1:indexLength,:)];
    
lineEnd =  [charAvg(end-indexLength:end,:);...
    interp1([1-aLenInterval,1],charAvg(end-1:end,1), 1.50,'linear','extrap'),...
    interp1([1-aLenInterval,1],charAvg(end-1:end,2), 1.50,'linear','extrap')];

%Find intercepts to divide line using Poly
[~,~,iIntStart] = polyxpoly(envelope(:,1),envelope(:,2),...
    lineStart(:,1),lineStart(:,2));
iIntStart = iIntStart(1);

[~,~,iIntEnd] = polyxpoly(envelope(:,1),envelope(:,2),...
    lineEnd(:,1),lineEnd(:,2));
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

% Resample corridors. Use nResamplePoints. Because corridors are
% non-monotonic, arc-length method discussed above is used. 
% Start with inner corridor. Magnitudes are being normalized.
segments = sqrt(((innerCorr(1:end-1,1)-innerCorr(2:end,1))./max(innerCorr(:,1))).^2 ...
    + ((innerCorr(1:end-1,2)-innerCorr(2:end,2))./max(innerCorr(:,2))).^2);
alen = cumsum([0;segments]);
alenResamp = linspace(0,max(alen),nvArg.nResamplePoints)';
innerCorr = [interp1(alen,innerCorr(:,1),alenResamp),...
    interp1(alen,innerCorr(:,2),alenResamp)];
% Outer Corridor
segments = sqrt(((outerCorr(1:end-1,1)-outerCorr(2:end,1))./max(outerCorr(:,1))).^2 ...
    + ((outerCorr(1:end-1,2)-outerCorr(2:end,2))./max(outerCorr(:,2))).^2);
alen = cumsum([0;segments]);
alenResamp = linspace(0,max(alen),nvArg.nResamplePoints)';
outerCorr = [interp1(alen,outerCorr(:,1),alenResamp),...
    interp1(alen,outerCorr(:,2),alenResamp)];

%% Draw Ellipses for debug
if strcmp(nvArg.Diagnostics,'on')
    % Plot corridors, avgs
    scatter(xx(:),yy(:),12,zz(:)>=1,'filled')
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
   