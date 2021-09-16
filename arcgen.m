%% ARCGen - Arc-length Response Corridor Generator
% 
% Created By:     D.C. Hartlen, M.ASc, EIT
% Date:           27-Jun-2021
% Updated By:     D.C. Hartlen, M.ASc, EIT
% Date:           19-Sep-2021
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
    arcgen(responseCurves,varargin)

%% Setup Name-Value Argument parser
nvArgObj = inputParser;
addParameter(nvArgObj, 'nResamplePoints',   100);
addParameter(nvArgObj, 'Diagnostics',       'off');
addParameter(nvArgObj, 'InvalidCurves',     []);
addParameter(nvArgObj, 'CorridorScaleFact', 1);
addParameter(nvArgObj, 'NormalizeCurves',   'on');
addParameter(nvArgObj, 'DeviationFact',     2);
addParameter(nvArgObj, 'EllipseKFact',      1);
addParameter(nvArgObj, 'CorridorRes',       100);
addParameter(nvArgObj, 'MinCorridorWidth',  0); 
addParameter(nvArgObj, 'nWarpCtrlPts',      0);
addParameter(nvArgObj, 'WarpingPenalty',    1e-2);
nvArgObj.KeepUnmatched = true;
nvArgObj.CaseSensitive = false;

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

% Error handling if NormalizeCurves argument is not defined correctly
else
    error('Normalization method not  recognized')
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

%% Resample response curve based on normalized arc-length
for iCurve=1:length(responseCurves)
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

for iPoints=1:nvArg.nResamplePoints
    clear temp; % probably cleaner way to do this.
    % collect specific point from each data curve
    for iCurve=1:length(responseCurves)
        temp(iCurve,:) = responseCurves(iCurve).normalizedCurve(iPoints,2:3);
    end
    charAvg(iPoints,:) = mean(temp,1);
    stdevData(iPoints,:) = std(temp,1);
end


%% Align normalized arc-length signals based on minimized correlation. 
% Enabled by option 'nWarpCtrlPts'. If 0, skip alignment.
if nvArg.nWarpCtrlPts > 0
    
    % Assemble signal matrices prior to correlation
    signalX = zeros(nvArg.nResamplePoints, length(responseCurves));
    signalY = zeros(nvArg.nResamplePoints, length(responseCurves));
    for i=1:length(responseCurves)
        signalX(:,i) = responseCurves(i).normalizedCurve(:,2);
        signalY(:,i) = responseCurves(i).normalizedCurve(:,3);
    end
    [meanCorrScore, corrArray] = evalCorrScore(signalX,signalY);
    
    % Optimize warp points for arbitrary n warping points. Build bounds,
    % constraints, and x0s
    nWarp = nvArg.nWarpCtrlPts;
    nSignal = length(responseCurves);
    
    if nWarp == 1   % nWarp == 1 is a special case as inequalites aren't needed
        x0 = 0.50.*ones(nSignal*2,1);
        lb = 0.15.*ones(nSignal*2,1);
        ub = 0.85.*ones(nSignal*2,1);
        A = [];
        b = [];
    elseif nWarp >= 15
        error('Specifying more than 10 interior warping points is not supported')
    else
        x0 = zeros(nWarp*(nSignal*2),1);
        for i = 1:nWarp
            x0(((i-1)*nSignal)+(1:nSignal) + (i-1)*nSignal) = i/(nWarp+1).*ones(nSignal,1);
            x0(((i-1)*nSignal)+(1:nSignal)+i*(nSignal)) = i/(nWarp+1).*ones(nSignal,1);
        end
        lb = 0.05.*ones(nWarp*(nSignal*2),1);
        ub = 0.95.*ones(nWarp*(nSignal*2),1);
        A = zeros((nWarp-1)*(nSignal*2), nWarp*(nSignal*2));
        b = -0.05.*ones((nWarp-1)*(nSignal*2), 1); % Force some separation between warped points
        for iSignal = 1:(nSignal*2)
            for iWarp = 1:(nWarp-1)
                A(iSignal+(iWarp-1)*(nSignal*2), iSignal+(iWarp-1)*(nSignal*2)) = 1;
                A(iSignal+(iWarp-1)*(nSignal*2), iSignal+iWarp*(nSignal*2)) = -1;
            end
        end
    end
    
    % Setup optimization options
    optOptions = optimoptions('fmincon',...
        'MaxFunctionEvaluations',max(3000, (nWarp+1).*1000),...
        'Display','off');
    
    % Execute optimization and compute warped signals
    optWarpArray = fmincon(@(x)warpingObjective(x,nWarp,...
        responseCurves,nvArg),...
        x0, A, b, [], [], lb, ub, [], optOptions);
    optWarpArray = reshape(optWarpArray,[],nWarp);
    [warpedSignals, signalX, signalY] = ...
        warpArcLength(optWarpArray,responseCurves,nvArg);
    

    % Compute correlation score
    [meanCorrScore, corrArray] = evalCorrScore(signalX,signalY);
    
    % Replace 'normalizedCurve' in 'responseCurve' and compute average and
    % standard deviation.
    for iSignal = 1:length(responseCurves)
        responseCurves(iSignal).normalizedCurve = warpedSignals{iSignal};
        responseCurves(iSignal).warpControlPoints = ...
            [[0,optWarpArray(iSignal+nSignal,:),1];...
            [0,optWarpArray(iSignal,:),1]];
    end
    for iPoints=1:nvArg.nResamplePoints
        clear temp; % probably cleaner way to do this.
        % collect specific point from each data curve
        for iSignal=1:length(responseCurves)
            temp(iSignal,:) = ...
                responseCurves(iSignal).normalizedCurve(iPoints,2:3);
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
if strcmp(nvArg.Diagnostics,'on') || strcmp(nvArg.Diagnostics,'detailed')
    figure('Name','Diagnostic Curves');
    cmap = lines(length(responseCurves));
    
    % Plot normalized x,y data
    subplot(2,2,1); hold on;
    for iCurve=1:length(responseCurves)
        pCurve(iCurve) = plot(responseCurves(iCurve).normalizedCurve(:,2),...
            responseCurves(iCurve).normalizedCurve(:,3),'.-',...
            'color',cmap(iCurve,:),...
            'DisplayName',responseCurves(iCurve).specId);
        if (strcmp(nvArg.NormalizeCurves,'off') || ...
                strcmp(nvArg.NormalizeCurves,'on'))
            continue
        else
            plot(responseCurves(iCurve).data(responseCurves(iCurve).alignInd,1),...
            responseCurves(iCurve).data(responseCurves(iCurve).alignInd,2),...
            'kx','LineWidth',2.0)
        end
    end
    xlabel('x-data')
    ylabel('y-data')
    legend(pCurve, 'location', 'Best')
    title('Arc-length Discretized Normalized Curves')
    
    % Plot warpping functions
    subplot(2,2,2); hold on
    clear pCurve
    if nvArg.nWarpCtrlPts > 0
        colours = lines(nSignal);
        for iSignal = 1:nSignal
            pCurve(iSignal) = plot(responseCurves(iSignal).data(:,4),...
                pchip([0,optWarpArray(iSignal+nSignal,:),1],[0,optWarpArray(iSignal,:),1],...
                responseCurves(iSignal).data(:,4)),...
                '.-','DisplayName',responseCurves(iSignal).specId,...
                'color',colours(iSignal,:),...
                'DisplayName',responseCurves(iCurve).specId);
            plot([0,optWarpArray(iSignal+nSignal,:),1],[0,optWarpArray(iSignal,:),1],'x',...
                'color',colours(iSignal,:),'MarkerSize',12,'LineWidth',2.0)
            title('Warping functions');
            legend(pCurve, 'location', 'Best')
        end
    else
        title('No Warping Performed');
    end
    plot([0,1],[0,1],'--','color',0.3.*[1,1,1])
    xlabel('Unwarped Normalized Arc-length')
    ylabel('Warped Normalized Arc-length')    
    
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

if strcmp(nvArg.Diagnostics,'detailed')
    % Plot ellipses    
    figure('Name','Ellipses and Corridor Extraction Debug'); hold on;
    % Scatter plot for debug
    cmap = cbrewer2('set2',2);
    colormap(cmap);
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
        zz(iPt,jPt) = max(...
            (((xx(iPt,jPt) - charAvg(:,1)).^2 ./ ...
            (stdevData(:,1).*nvArg.EllipseKFact).^2 ...
            + (yy(iPt,jPt) - charAvg(:,2)).^2 ./ ...
            (stdevData(:,2).*nvArg.EllipseKFact).^2).^-1));
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
% Initalize line segments for speed. This line may cause issues, as it
% assumes maximum size. Bump up 10 if it does. 
lineSegments = zeros(10*max(nvArg.nResamplePoints,nvArg.CorridorRes),4); 
iSeg = 0;
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
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))];        
            case 3
                % West-North
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)];
            case 4
                % North-South
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt) ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1)),yy(iPt+1,jPt)];
            case 5
                % North-East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))];        
            case 6  % Ambiguous 
                centerVal = mean([zz(iPt,jPt), zz(iPt+1,jPt), zz(iPt+1,jPt+1), zz(iPt, jPt+1)]);
                if centerVal >= 1
                    % West-North
                    iSeg = iSeg+1;
                    lineSegments(iSeg,:) = ...
                        [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                        interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)];
                    % South - East
                    iSeg = iSeg+1;
                    lineSegments(iSeg,:) = ...
                        [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                        xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))];
                else
                    % South-West
                    iSeg = iSeg+1;
                    lineSegments(iSeg,:) = ...
                        [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                        xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))];
                    % North-East
                    iSeg = iSeg+1;
                    lineSegments(iSeg,:) = ...
                        [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                        xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))];
                end
            case 7
                % West-East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [xx(iPt,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))];
            case 8
                % South - East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))];  
            case 9
                % South - East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))];
            case 10
                % West-East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [xx(iPt,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))];
            case 11 % Ambiguous
                centerVal = mean([zz(iPt,jPt), zz(iPt+1,jPt), zz(iPt+1,jPt+1), zz(iPt, jPt+1)]);
                if centerVal >= 1
                    % South-West
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                        [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                        xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))];
                    % North-East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                        [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                        xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))];
                else
                    % West-North
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                        [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                        interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)];
                    % South-East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                        [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                        xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))];
                end
            case 12
                % North-East
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))];
            case 13
                % North-South
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt) ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1)),yy(iPt+1,jPt)];
            case 14
                % West-North
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)];
            case 15
                % South-West
                iSeg = iSeg+1;
                lineSegments(iSeg,:) = ...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]; 
            case 16
                % No vertices 
        end
    end
end
lineSegments = lineSegments(1:iSeg,:);

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
if strcmp(nvArg.Diagnostics,'detailed')
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

aLenExtension = abs(aLenInterval./(charAvg(1,:)-charAvg(2,:)))...
    .*1.1.*max(stdevData);
aLenExtension(isinf(aLenExtension)) = 0;
aLenExtension = max(aLenExtension);

lineStart = [...
    interp1([0,aLenInterval],charAvg(1:2,1), -aLenExtension,'linear','extrap'),...
    interp1([0,aLenInterval],charAvg(1:2,2), -aLenExtension,'linear','extrap');...
    charAvg(1:indexLength,:)];

lineEnd =  [charAvg(end-indexLength:end,:);...
    interp1([1,1-aLenInterval],[charAvg(end,1),charAvg(end-1,1)],...
    (1+aLenExtension),'linear','extrap'),...
    interp1([1,1-aLenInterval],[charAvg(end,2),charAvg(end-1,2)],...
    (1+aLenExtension),'linear','extrap')];

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

%% Draw extension lines and sampling points to MS plot
if strcmp(nvArg.Diagnostics,'detailed')
    % Plot corridors, avgs
    scatter(xx(:),yy(:),12,zz(:)>=1,'filled')
    plot(lineStart(:,1),lineStart(:,2),'.-k','DisplayName','Char Avg',...
        'LineWidth',2.0,'MarkerSize',16)
    plot(lineEnd(:,1),lineEnd(:,2),'.-k','DisplayName','Char Avg',...
        'LineWidth',2.0,'MarkerSize',16)
    xlim([min(xx(:)),max(xx(:))])
    ylim([min(yy(:)),max(yy(:))])
end
processedCurveData = responseCurves;
end  % End main function

%% helper function to perform linear interpolation to an isovalue of 1 only
function val = interpVal(x1, y1, x2, y2)
    val = x1+(x2-x1)*(1-y1)/(y2-y1);
end
   
%% Function used to evaluate correlation score between signals
function [meanCorrScore, corrScoreArray] = evalCorrScore(signalsX,signalsY)
% Correlation score taken from the work of Nusholtz et al. (2009)
% Compute cross-correlation matrix of all signals to each other
corrMatX = corrcoef(signalsX);
corrMatY = corrcoef(signalsY);
% Convert matrices to a single score
nCurve = size(corrMatX,2);
corrScoreX = (1/(nCurve*(nCurve-1)))*(sum(sum(corrMatX))-nCurve);
corrScoreY = (1/(nCurve*(nCurve-1)))*(sum(sum(corrMatY))-nCurve);
% Compute a single metric for optimization purposes. Using simple mean
meanCorrScore = 0.5*(corrScoreX+corrScoreY);
corrScoreArray = [corrScoreX, corrScoreY];
end

%% Function used to compute objective for optimization
function [optScore, penaltyScore] = warpingObjective(optimWarp,nCtrlPts,responseCurves,nvArg)
% Control points are equally spaced in arc-length. 
% optimwarp is a column vector with first warped control point in the
% first nCurve indices, then 2nd control point in the next nCurve indices

% warpArray = reshape(optimWarp,length(responseCurves),nCtrlPts);
nSignal = length(responseCurves);
warpArray = reshape(optimWarp,[],nCtrlPts);
% Compute a warping penalty
penaltyScore = warpingPenalty(warpArray,nvArg.WarpingPenalty,nvArg);
penaltyScore = mean(penaltyScore);

% Perform warping
[~, signalsX, signalsY] = warpArcLength(warpArray,responseCurves,nvArg);
% Compute correlation score
[corrScore, ~] = evalCorrScore(signalsX,signalsY);
% corrScore is a maximization goal. Turn into a minimization goal
optScore = 1-corrScore+penaltyScore;

end

%% Function used to warp arc-length
function [warpedSignals, signalsX, signalsY]...
    = warpArcLength(warpArray, responseCurves, nvArg)
% Warp array: each row is warping points for an input signal, each column
% is warped point. Control points are interpolated  on [0,1] assuming
% equal spacing. 
[~, nCtrlPts] = size(warpArray);
nCurves = length(responseCurves);


% lmCtrlPts = linspace(0,1,2+nCtrlPts);
% lmCtrlPts = [0,warpArray(end,:),1];

% Initialize matrices
signalsX = zeros(nvArg.nResamplePoints, nCurves);
signalsY = zeros(nvArg.nResamplePoints, nCurves);
warpedSignals = cell(nCurves,1);

for iCurve = 1:nCurves
    % Assign responseCurve data array to matrix for brevity
    curve = responseCurves(iCurve).data;
    
    lmCtrlPts = [0,warpArray(iCurve+nCurves,:),1];
    
    % prepend 0 and append 1 to warp points for this curve to create valid
    % control points. 
    warpedCtrlPts = [0,warpArray(iCurve,:),1];
    
    % Construct warping function using SLM. This warps lmAlen to shiftAlen.
    % Use warping fuction to map computed arc-lengths onto the shifted
    % system. use built-in pchip function. This is a peicewise monotonic 
    % cubic spline. Signifincantly faster than SLM. 
    warpedNormAlen = pchip(lmCtrlPts,warpedCtrlPts,curve(:,4));
      
    % Now uniformly resample normalzied arc-length
    resamNormwarpedAlen = linspace(0,1, nvArg.nResamplePoints)';
    resampX = interp1(warpedNormAlen, curve(:,1), resamNormwarpedAlen,'linear','extrap');
    resampY = interp1(warpedNormAlen, curve(:,2), resamNormwarpedAlen,'linear','extrap');
    % Assign to array for correlation calc
    signalsX(:,iCurve) = resampX;
    signalsY(:,iCurve) = resampY;
    
    % Assemble a cell array containing arrays of resampled signals. Similar
    % to 'normalizedCurve' in 'responseCurves' structure
    warpedSignals{iCurve} = [resamNormwarpedAlen,resampX,resampY];
end

end

%% Penalty function to prevent plateaus and extreme divergence in warping functions
function [penaltyScores] = warpingPenalty(warpArray,penaltyFactor,nvArg)
% Compute an array of penalty scores based on MSE between linear, unwarped
% arc-length and warped arc-length. Aim is to help prevent plateauing. 
[nCurves, nCtrlPts] = size(warpArray);
nCurves = nCurves/2;
% lmCtrlPts = [0, warpArray(end,:), 1];
penaltyScores = zeros(nCurves,1);
unwarpedAlen = linspace(0,1,nvArg.nResamplePoints);

for iCurve=1:nCurves
    penaltyScores(iCurve) = sum((unwarpedAlen - ...
        pchip([0,warpArray(iCurve+nCurves,:),1],...
        [0,warpArray(iCurve,:),1],unwarpedAlen)).^2);
end

penaltyScores = penaltyScores.*penaltyFactor;
end
    
