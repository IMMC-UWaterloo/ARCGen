%% ARCGen - Arc-length Response Corridor Generator
%
% ARCGen, short for Arc-length Response Corridor Generation, provides
% automated calculation of a characteristic average and response corridors
% on input signals regardless of if said signals are non-monotonic or
% hystertic. This is accomplished by re-parameterizing input signals based
% on arc-length. Corridors are extracted using a marching squares
% algorithm.
%
% If you use ARCGen in your research, please use the following citation:
%     Hartlen D.C. and Cronin D.S. (2022), "Arc-Length Re-Parametrization 
%        and Signal Registration to Determine a Characteristic Average and
%        Statistical Response Corridors of Biomechanical Data." Frontiers 
%        in Bioengineering and Biotechnology* 10:843148. 
%        doi: 10.3389/fbioe.2022.843148
%
% ARCGen is released under a GNU GPL v3 license. No warranty or support is
% provided. The authors hold no responsibility for the validity, accuracy, 
% or applicability of any results obtained from this code.
%
% This function has one mandatory input, four outputs, and many optional
% inputs. Optional inputs are defined using name-value pair arguments. 
%
% Usage notes: 
% It is common to see errors when running this function if the number of
% resampling points or corridor resolution is too sparse or signals 
% exhibit significant variablity not accounted for through signal 
% registration. This tends to manifest in either truncated corridors or the
% code termininating in an error. Often increasing resampling points or
% corridor resolution. Turning 'Diagnostics' to 'detailed' can help 
% identify these issues. 
%
% Computed corridors will often not extend all the way to the shared origin
% of input signals. This is because small low st. dev. at this shared point
% is too low to be captured during corridor extraction with the marching
% squares algorithm. There are two solutions to this problem. First, one
% could force a minimum corridors size using the 'MinCorridorWidth' option.
% Second, one could manually extend corridors in post-processing. 
%
% MANDATORY INPUTS:
% -----------------
% inputSignals: ARCGen can accomadate three types of input format
% 1) A [nSignal,2] structured array consisting of the following
%       entries. Entries are case-senstive
%   + data: an [m,2] array containing ordered x-y data
%   + specId: character array containing an identifier for each signal
% 2) A [nSignal,1] structured array consisting of only signal data, no 
%       signal IDs. Entries are case-senstive
%   + data: an [m,2] array containing ordered x-y data
% 3) A cell array of length nSignal containing [m,2] arrays of each input
%       signal. 
% 
% OPTIONAL INPUTS:
% ----------------
% nResamplePoints: integer defining the number of points used to
%       re-parameterize input signals. Default: 100. 
% CorridorRes: integeer defining the number of grid points used for the
%       marching squares algorithm. The sampling grid for the marching
%       squares algorithm extends 120% of extreme corridors. This parameter
%       defines the number of points along each side of the grid. 
%       Default: 100. It is common to increase this significantly. 
% NormalizeSignals: character arry used to turn on signal normalization.
%       Options: 'on' (default), 'off'
% EllipseKFact: float used to scale the major and minor axis of the
%       ellipses used for corridor generation. This value corrisponds to 
%       the square root of the chi-squared CDF. Default: 1.0 (creates 
%       corridors one standard deviation along the x and y axes)
% Diagnostics: character array used to activate diagnostic plots. Useful
%       for debugging errors. Options: 'off' (default), 'on', 'detailed'. 
% MinCorridorWidth: Factor used to enforce a minimum corridor width. Any
%       st.dev. less than 'MinCorridorFactor'*max(st.dev.) is replaced with
%       'MinCorridorFactor'*max(st.dev.). x & y axes are handled
%       separately. A value of 0 (default) disables forcing minimum width.
% nWarpCtrlPts: integer that sets the number of interior control points
%       used for signal registration. A value of 0 (default) disables
%       signal registration
% WarpingPenalty: float specifying the penalty factor used during the
%       signal registration process. A value of 10^-2 (default) to 10^3 is
%       recommended, but the exact value will need to be tuned to a
%       specific problem. 
% UseParrallel: Character array used to enable parallel thread calculations
%       for signal registration and envelope extraction. Significantly
%       reduces runtime when signals have 100k+ points or 500+ resampling
%       points and corridor resolution. Requires the Parallel Computing 
%       Toolbox Options: 'on', 'off' (default). 
%
% MANDATORY OUTPUTS:
% ------------------
% charAvg: an [nResamplePoints,2] array containing the computed
%       characteristic average.
% innerCorr: an [nResamplePoints,2] array containing points defining the
%       inner corridor
% outerCorr: an [nResamplePoints,2] array containing points defining the 
%       outer corridor
%
% OPTIONAL OUTPUTS:
% ----------------
% processedSignalData: a structure array that outputs processed signals,
%       basic statistics, and warping control poitns
% debugData: a structure that provides a wealth of debugging information,
%       including raw average and st. dev. data, correlation scores before
%       and after registration, and other. 
%
% Copyright (c) 2022 Devon C. Hartlen

function [charAvg, innerCorr, outerCorr, varargout] = ...
    arcgen(inputSignals,varargin)

%% Setup Name-Value Argument parser
nvArgObj = inputParser;
addParameter(nvArgObj, 'nResamplePoints',   100);
addParameter(nvArgObj, 'Diagnostics',       'off');
addParameter(nvArgObj, 'NormalizeSignals',  'on');
addParameter(nvArgObj, 'EllipseKFact',      1);
addParameter(nvArgObj, 'CorridorRes',       100);
addParameter(nvArgObj, 'MinCorridorWidth',  0); 
addParameter(nvArgObj, 'nWarpCtrlPts',      0);
addParameter(nvArgObj, 'WarpingPenalty',    1e-2);
addParameter(nvArgObj, 'UseParallel',       'off');
nvArgObj.KeepUnmatched = true;
nvArgObj.CaseSensitive = false;

parse(nvArgObj,varargin{:});
nvArg = nvArgObj.Results;  % Structure created for convenience

% check if parallel toobox is installed, then if parpool is running. Error
% out if not installed, start pool if not already started. 
v = ver;
hasParallel = any(strcmp(cellstr(char(v.Name)), 'Parallel Computing Toolbox'));
if strcmp(nvArg.UseParallel,'on')
    if ~hasParallel
        error('Parallel Computing Toolbox is not installed. Set option UseParallel to off')
    end
    p = gcp('nocreate');
    if isempty(p)
        parpool();
    end
end

%% Add third party functions to path
funcPath = mfilename('fullpath');
funcPath = fileparts(funcPath);
addpath(fullfile(funcPath, 'ThirdPartyFunctions'));

%% Process input options
% Check if structure with specID, struct w/o specID, cell array. Error out
% otherwise. Places inputs into structure format.
if isstruct(inputSignals)
    if ~isfield(inputSignals,'specId')
        for iSignal = 1:length(inputSignals)
            inputSignals(iSignal).specId = ...
                ['Signal ' num2str(iSignal,'%3d')];
        end
    end
elseif iscell(inputSignals)
    inputSignals = cell2struct(inputSignals,'data');
    for iSignal = 1:length(inputSignals)
       inputSignals(iSignal).specId = ['Signal ' num2str(iSignal,'%3d')];
    end
end

%% Compute arclength based on input signal datapoints
% Do not perform normalization
if strcmp(nvArg.NormalizeSignals,'off')
    for iSignal = 1:length(inputSignals)
        temp = inputSignals(iSignal).data; % Temporary for conveinence
        % Compute arc-length between each data point
        segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
            + (temp(1:end-1,2)-temp(2:end,2)).^2);
        alen = cumsum([0;segments]);
        % Append cumulative arc length to data array
        inputSignals(iSignal).data = [inputSignals(iSignal).data,alen];
        % Compute normalized arc-length
        inputSignals(iSignal).maxAlen = max(alen);
        inputSignals(iSignal).data = [inputSignals(iSignal).data,...
            alen./inputSignals(iSignal).maxAlen];
        % Determine max [x,y] data
        tempMax = max(temp,[],1);
        inputSignals(iSignal).xMax = tempMax(1);
        inputSignals(iSignal).yMax = tempMax(2);
        % Remove spurious duplicates
        [~,index,~] = unique(inputSignals(iSignal).data(:,4));
        inputSignals(iSignal).data = inputSignals(iSignal).data(index,:);
    end
    
% Perform magnitude normalization based on bounding box
elseif strcmp(nvArg.NormalizeSignals,'on')
    % Determine bounding box of individual signals
    for iSignal = 1:length(inputSignals)
        tempMin = min(inputSignals(iSignal).data,[],1);
        inputSignals(iSignal).xMin = tempMin(1);
        inputSignals(iSignal).yMin = tempMin(2);
        tempMax = max(inputSignals(iSignal).data,[],1);
        inputSignals(iSignal).xMax = tempMax(1);
        inputSignals(iSignal).yMax = tempMax(2);
    end
    xBound = [mean([inputSignals.xMin]), mean([inputSignals.xMax])];
    yBound = [mean([inputSignals.yMin]), mean([inputSignals.yMax])];
    % Normalize the axis of each signal, then do arc-length calcs
    for iSignal = 1:length(inputSignals)
        temp = inputSignals(iSignal).data; % Temporary for conveinence
        % Normalize from bounding box to [-1,1]
        temp = [temp(:,1)./(xBound(2)-xBound(1)),...
            temp(:,2)./(yBound(2)-yBound(1))];
        % Compute arc-length between each data point
        segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
            + (temp(1:end-1,2)-temp(2:end,2)).^2);
        alen = cumsum([0;segments]);
        % Append cumulative arc length to data array
        inputSignals(iSignal).data = [inputSignals(iSignal).data,alen];
        % Compute normalized arc-length
        inputSignals(iSignal).maxAlen = max(alen);
        inputSignals(iSignal).data = [inputSignals(iSignal).data,...
            alen./inputSignals(iSignal).maxAlen];
        % Determine max [x,y] data
        tempMax = max(abs(temp),[],1);
        inputSignals(iSignal).xNormMax = tempMax(1);
        inputSignals(iSignal).yNormMax = tempMax(2);
        % Remove spurious duplicates
        [~,index,~] = unique(inputSignals(iSignal).data(:,4));
        inputSignals(iSignal).data = inputSignals(iSignal).data(index,:);
    end

% Error handling if NormalizeSignals argument is not defined correctly
else
    error('Normalization method not  recognized')
end

% Compute mean and median arc-length deviation
meanAlen = mean([inputSignals.maxAlen]);
for iSignal=1:length(inputSignals)
    inputSignals(iSignal).meanDevs = ...
        inputSignals(iSignal).maxAlen-meanAlen;
end

medianAlen = median([inputSignals.maxAlen]);
for iSignal=1:length(inputSignals)
    inputSignals(iSignal).medianDev = ...
        inputSignals(iSignal).maxAlen-medianAlen;
end

%% Resample response signal based on normalized arc-length
for iSignal=1:length(inputSignals)
    % Linear-interpolation for x,y data against arc-length
    normAlen = linspace(0,inputSignals(iSignal).data(end,4),...
        nvArg.nResamplePoints)';
    resampX = interp1(inputSignals(iSignal).data(:,4),...
        inputSignals(iSignal).data(:,1), normAlen);
    resampY = interp1(inputSignals(iSignal).data(:,4),...
        inputSignals(iSignal).data(:,2), normAlen);
    % Resulting array is normalized arc-length, resampled x, resam. y
    inputSignals(iSignal).normalizedSignal = [normAlen, resampX, resampY];
end
    
%% For each resampled point, determine average and standard deviation across signals
% Initialize arrays
charAvg = zeros(nvArg.nResamplePoints,2);
stdevData = zeros(nvArg.nResamplePoints,2);

for iPoints=1:nvArg.nResamplePoints
    clear temp; % probably cleaner way to do this.
    % collect specific point from each signal
    for iSignal=1:length(inputSignals)
        temp(iSignal,:) = inputSignals(iSignal).normalizedSignal(iPoints,2:3);
    end
    charAvg(iPoints,:) = mean(temp,1);
    stdevData(iPoints,:) = std(temp,1);
end
% Assign characteristic average and st. dev. data to a debug structure
debugOutput.charAvg = charAvg;
debugOutput.stdevData = stdevData;

%% Align normalized arc-length signals based on minimized correlation. 
% Enabled by option 'nWarpCtrlPts'. If 0, skip alignment.
if nvArg.nWarpCtrlPts > 0
    % Assemble signal matrices prior to correlation
    signalX = zeros(nvArg.nResamplePoints, length(inputSignals));
    signalY = zeros(nvArg.nResamplePoints, length(inputSignals));
    for i=1:length(inputSignals)
        signalX(:,i) = inputSignals(i).normalizedSignal(:,2);
        signalY(:,i) = inputSignals(i).normalizedSignal(:,3);
    end
    [meanCorrScore, corrArray] = evalCorrScore(signalX,signalY);
    % Assign pre-optimized correlation scores to debug structure
    debugOutput.preWarpCorrArray = corrArray;
    debugOutput.preWarpMeanCorrScore = meanCorrScore;
    
    % Optimize warp points for arbitrary n warping points. Build bounds,
    % constraints, and x0s
    nWarp = nvArg.nWarpCtrlPts;
    nSignal = length(inputSignals);
    
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
    
    % Setup optimization options ('UseParallel' option active here)
    if strcmp(nvArg.UseParallel,'on')
        optOptions = optimoptions('fmincon',...
            'MaxFunctionEvaluations',max(3000, (nWarp+1).*1000),...
            'Display','off',...
            'UseParallel',true);
    else
        optOptions = optimoptions('fmincon',...
        'MaxFunctionEvaluations',max(3000, (nWarp+1).*1000),...
        'Display','off',...
        'UseParallel',false);
    end
    
    % Execute optimization and compute warped signals
    optWarpArray = fmincon(@(x)warpingObjective(x,nWarp,...
        inputSignals,nvArg),...
        x0, A, b, [], [], lb, ub, [], optOptions);
    optWarpArray = reshape(optWarpArray,[],nWarp);
    [warpedSignals, signalX, signalY] = ...
        warpArcLength(optWarpArray,inputSignals,nvArg.nResamplePoints);
    

    % Compute correlation score
    [meanCorrScore, corrArray] = evalCorrScore(signalX,signalY);
    % Assign warped correlation scores to debug structure
    debugOutput.warpedCorrArray = corrArray;
    debugOutput.warpedMeanCorrScore = meanCorrScore;
    
    % Replace 'normalizedSignal' in 'responseSignal' and compute average and
    % standard deviation.
    for iSignal = 1:length(inputSignals)
        inputSignals(iSignal).normalizedSignal = warpedSignals{iSignal};
        inputSignals(iSignal).warpControlPoints = ...
            [[0,optWarpArray(iSignal+nSignal,:),1];...
            [0,optWarpArray(iSignal,:),1]];
    end
    for iPoints=1:nvArg.nResamplePoints
        clear temp; % probably cleaner way to do this.
        % collect specific point from each signal
        for iSignal=1:length(inputSignals)
            temp(iSignal,:) = ...
                inputSignals(iSignal).normalizedSignal(iPoints,2:3);
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

%% Diagnostic: Plot normalized signals and St. Devs. 
if strcmp(nvArg.Diagnostics,'on') || strcmpi(nvArg.Diagnostics,'detailed')
    figure('Name','Diagnostic Signals');
    cmap = lines(length(inputSignals));
    
    % Plot normalized x,y data
    subplot(2,2,1); hold on;
    for iSignal=1:length(inputSignals)
        pSignal(iSignal) = plot(inputSignals(iSignal).normalizedSignal(:,2),...
            inputSignals(iSignal).normalizedSignal(:,3),'.-',...
            'color',cmap(iSignal,:),...
            'DisplayName',inputSignals(iSignal).specId);
        if (strcmp(nvArg.NormalizeSignals,'off') || ...
                strcmp(nvArg.NormalizeSignals,'on'))
            continue
        else
            plot(inputSignals(iSignal).data(inputSignals(iSignal).alignInd,1),...
            inputSignals(iSignal).data(inputSignals(iSignal).alignInd,2),...
            'kx','LineWidth',2.0)
        end
    end
    xlabel('x-data')
    ylabel('y-data')
    legend(pSignal, 'location', 'Best')
    title('Arc-length Discretized Normalized Signals')
    
    % Plot warpping functions
    subplot(2,2,2); hold on
    clear pSignal
    if nvArg.nWarpCtrlPts > 0
        colours = lines(nSignal);
        for iSignal = 1:nSignal
            pSignal(iSignal) = plot(inputSignals(iSignal).data(:,4),...
                pchip([0,optWarpArray(iSignal+nSignal,:),1],[0,optWarpArray(iSignal,:),1],...
                inputSignals(iSignal).data(:,4)),...
                '.-','DisplayName',inputSignals(iSignal).specId,...
                'color',colours(iSignal,:),...
                'DisplayName',inputSignals(iSignal).specId);
            plot([0,optWarpArray(iSignal+nSignal,:),1],[0,optWarpArray(iSignal,:),1],'x',...
                'color',colours(iSignal,:),'MarkerSize',12,'LineWidth',2.0)
            title('Warping functions');
            legend(pSignal, 'location', 'Best')
        end
    else
        title('No Warping Performed');
    end
    plot([0,1],[0,1],'--','color',0.3.*[1,1,1])
    xlabel('Unwarped Normalized Arc-length')
    ylabel('Warped Normalized Arc-length')    
    
    % Plot normalized x data against arc-length with st. dev.
    subplot(2,2,3); hold on;
    errorbar(inputSignals(1).normalizedSignal(:,1),charAvg(:,1),...
        stdevData(:,1),'color',0.5.*[1,1,1])
    cmap = lines;
    for iSignal=1:length(inputSignals)
        plot(inputSignals(iSignal).normalizedSignal(:,1),...
            inputSignals(iSignal).normalizedSignal(:,2),'.-',...
            'color',cmap(iSignal,:))
    end
    xlabel('Normalized Arc-length')
    ylabel('x-data')
    title('Average and St.Dev. of X-Data')
    
    % Plot normalized y data against arc-length with st. dev.
    subplot(2,2,4); hold on;
    errorbar(inputSignals(1).normalizedSignal(:,1),charAvg(:,2),...
        stdevData(:,2),'color',0.5.*[1,1,1])
    cmap = lines;
    for iSignal=1:length(inputSignals)
        plot(inputSignals(iSignal).normalizedSignal(:,1),...
            inputSignals(iSignal).normalizedSignal(:,3),'.-',...
            'color',cmap(iSignal,:))
    end
    xlabel('Normalized Arc-length')
    ylabel('y-data')
    title('Average and St.Dev. of Y-Data')
end

if strcmpi(nvArg.Diagnostics,'detailed')
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
    cmap = lines(length(inputSignals));
    for iSignal=1:length(inputSignals)
        plot(inputSignals(iSignal).data(:,1),...
            inputSignals(iSignal).data(:,2),'-',...
            'DisplayName',inputSignals(iSignal).specId,...
            'Color', cmap(iSignal,:))
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
kFact = nvArg.EllipseKFact; % faster if no struct call in inner loop. 
nRes = nvArg.CorridorRes;   % again, for speed
% If 'UseParallel' is 'on', grid evaluation is performed using a parallel
% for loop. 
if strcmp(nvArg.UseParallel,'on')
    parfor iPt = 1:nRes
        for jPt = 1:nRes
            zz(iPt,jPt) = max(...
                (((xx(iPt,jPt) - charAvg(:,1)).^2 ./ ...
                (stdevData(:,1).*kFact).^2 ...
                + (yy(iPt,jPt) - charAvg(:,2)).^2 ./ ...
                (stdevData(:,2).*kFact).^2).^-1));
        end
    end
% otherwise, use a standard forloop
else
    for iPt = 1:nRes
        for jPt = 1:nRes
            zz(iPt,jPt) = max(...
                (((xx(iPt,jPt) - charAvg(:,1)).^2 ./ ...
                (stdevData(:,1).*kFact).^2 ...
                + (yy(iPt,jPt) - charAvg(:,2)).^2 ./ ...
                (stdevData(:,2).*kFact).^2).^-1));
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

% Extract list of unique vertices from line segmens
vertices = [lineSegments(:,1:2);lineSegments(:,3:4)];
vertices = uniquetol(vertices,eps,'ByRows', true);

% Create a vertex connectivity table. The 1e-12 value is here because
% floats don't round well and == will not work. 
vertConn = zeros(size(lineSegments,1),2);
for i = 1:length(vertConn)
    index = all(abs(lineSegments(:,1:2) - vertices(i,:)) < 1e-12,2);
    vertConn(index,1) = i;
    index = all(abs(lineSegments(:,3:4) - vertices(i,:)) < 1e-12,2);
    vertConn(index,2) = i;
end

%% Start line segments sorting and envelope extraction
nEnvelopes = 1;
allEnvelopes(1,1) = 1;     % First entry is always vertex 1
for i = 1:size(vertConn,1)-1
    % save vertex to find 
    vertToFind = vertConn(i,2);
    j = i+1; % helper index
    % Find connecting node
    foundShiftedInd =...
        find(any(vertConn(j:end,:) == vertToFind,2), 1, 'first');
    % If we have found an index 
    if ~isempty(foundShiftedInd)
        foundInd = foundShiftedInd + i;
        % swap found vert conn row with j row
        temp = vertConn(j,:);
        % Now, decide whether to flip found row. We want vertex 2 of 
        % previous line to be node 1 of the new line. 
        if (vertConn(foundInd,1) == vertToFind)
            vertConn(j,:) = vertConn(foundInd, [1,2]);
        else 
            vertConn(j,:) = vertConn(foundInd, [2,1]);
        end
        % Logic to prevent overwriting, if found row is next row. 
        if (foundInd ~= j)
            vertConn(foundInd,:) = temp;
        end
    % If we did not find an index, we either may have an open envelope or
    % envelope may be convex and loops back on itself. 
    else
        % Check to see if we can find the first vertex in envelope
        % appearing again (check for closure)
        vertToFind = vertConn(allEnvelopes(nEnvelopes,1));
        foundShiftedInd = ...
            find(any(vertConn(j:end,:) == vertToFind,2), 1, 'first');
        % If we do not find an index, it means this envelope is complete
        % and manifold
        if isempty(foundShiftedInd)
            % Assign indices to finish current envelope, initialize next
            allEnvelopes(nEnvelopes,2) = i;
            nEnvelopes = nEnvelopes + 1;
            allEnvelopes(nEnvelopes, 1) = j;
        else
            % This error should only occur if envelopes extend beyond
            % sampling grid, which they should not. 
            error('Literal Edge Case')
        end
        
    end
end
allEnvelopes(nEnvelopes,2) = j;

% Find largest envelope
[~,envInds] = max(allEnvelopes(:,2)-allEnvelopes(:,1));

% Convert indices in evelopes to array of (x,y)
envInds = allEnvelopes(envInds, :);
envelope = vertices(vertConn(envInds(1):envInds(2),1),:);

% For debugging, plot all envelopes
if strcmpi(nvArg.Diagnostics,'detailed')
    for iEnv = 1:nEnvelopes
        envInds = allEnvelopes(iEnv, :);
        plot(vertices(vertConn(envInds(1):envInds(2),1),1),...
            vertices(vertConn(envInds(1):envInds(2),1),2),...
            '.-b','LineWidth',1.0)
    end
end

%% Divide the envelope into corridors. 
% To break the largest envelop into inner and outer corridors, we need to
% account for several edge cases. First, we test to see if there are any
% intercepts of the characteristic average and the largest envelope. 
closedEnvelope = [envelope; envelope(1,:)];
[~,~,indexIntercept] = polyxpoly(closedEnvelope(:,1),closedEnvelope(:,2),...
    charAvg(:,1),charAvg(:,2));

% If we find two intercepts, then we have no problem
if size(indexIntercept,1) >=2
    iIntStart = indexIntercept(1,1);
    iIntEnd = indexIntercept(end,1);

% If we find only one intercept, we need to determine if the intercept is a
% the start or end of the envelope. Then we need to extend the opposite
% side of the characteristic average to intercept the envelope. 
elseif size(indexIntercept,1) == 1
    % Compute extension 
    aLenInterval = 1./nvArg.nResamplePoints;
    indexLength = round(0.2*length(charAvg));
    
    aLenExtension = abs(aLenInterval./(charAvg(1,:)-charAvg(2,:)))...
        .*1.1.*max(stdevData);
    aLenExtension(isinf(aLenExtension)) = 0;
    aLenExtension = max(aLenExtension);
    % If the single found point is inside the envelope, the found intercept
    % is at the end. Therefore extend the start
    if inpolygon(charAvg(indexIntercept(2),1),...
            charAvg(indexIntercept(2),2), envelope(:,1),envelope(:,2))
        
        iIntEnd = indexIntercept(1);
        lineStart = [...
            interp1([0,aLenInterval],charAvg(1:2,1), -aLenExtension,'linear','extrap'),...
            interp1([0,aLenInterval],charAvg(1:2,2), -aLenExtension,'linear','extrap');...
            charAvg(1:indexLength,:)];

        % Plot line extensions for envelope splitting
        if strcmpi(nvArg.Diagnostics, 'Detailed')
            plot(lineStart(1:2,1), lineStart(1:2,2),'k:',...
                'DisplayName','Start Ext.')
        end

        %Find intercepts to divide line using Poly
        [~,~,iIntStart] = polyxpoly(closedEnvelope(:,1),closedEnvelope(:,2),...
            lineStart(:,1),lineStart(:,2));
        iIntStart = iIntStart(1);
        % If the single found point is outside the envelope, the found
    % intercept is the start
    else
        iIntStart = indexIntercept(1);
        lineEnd =  [charAvg(end-indexLength:end,:);...
            interp1([1,1-aLenInterval],[charAvg(end,1),charAvg(end-1,1)],...
            (1+aLenExtension),'linear','extrap'),...
            interp1([1,1-aLenInterval],[charAvg(end,2),charAvg(end-1,2)],...
            (1+aLenExtension),'linear','extrap')];

        % Plot line extensions for envelope splitting
        if strcmpi(nvArg.Diagnostics, 'Detailed')
            plot(lineEnd(end-1:end,1), lineEnd(end-1:end,2), 'k:', ...
                'DisplayName', 'End Ext.')
        end

        %Find intercepts to divide line using Poly
        [~,~,iIntEnd] = polyxpoly(closedEnvelope(:,1),closedEnvelope(:,2),...
            lineEnd(:,1),lineEnd(:,2));
        iIntEnd = iIntEnd(1);
    end
    
% If we find no intercepts, we need to extend both sides of characteristic
% average to intercept the envelop.
else
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

    % Plot line extensions for envelope splitting
    if strcmpi(nvArg.Diagnostics, 'Detailed')
        plot(lineStart(1:2,1), lineStart(1:2,2),'k:',...
            'DisplayName','Start Ext.')
        plot(lineEnd(end-1:end,1), lineEnd(end-1:end,2), 'k:', ...
            'DisplayName', 'End Ext.')
    end
    
    %Find intercepts to divide line using Poly
    [~,~,iIntStart] = polyxpoly(closedEnvelope(:,1),closedEnvelope(:,2),...
        lineStart(:,1),lineStart(:,2));
    iIntStart = iIntStart(1);
    
    [~,~,iIntEnd] = polyxpoly(closedEnvelope(:,1),closedEnvelope(:,2),...
        lineEnd(:,1),lineEnd(:,2));
    iIntEnd = iIntEnd(1);
end

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
if strcmpi(nvArg.Diagnostics,'detailed')
    % Plot corridors, avgs
    scatter(xx(:),yy(:),12,zz(:)>=1,'filled')
%     plot(lineStart(:,1),lineStart(:,2),'.-k','DisplayName','Char Avg',...
%         'LineWidth',2.0,'MarkerSize',16)
%     plot(lineEnd(:,1),lineEnd(:,2),'.-k','DisplayName','Char Avg',...
%         'LineWidth',2.0,'MarkerSize',16)
    xlim([min(xx(:)),max(xx(:))])
    ylim([min(yy(:)),max(yy(:))])
end

varargout{1} = inputSignals;
varargout{2} = debugOutput;
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
nSignal = size(corrMatX,2);
corrScoreX = (1/(nSignal*(nSignal-1)))*(sum(sum(corrMatX))-nSignal);
corrScoreY = (1/(nSignal*(nSignal-1)))*(sum(sum(corrMatY))-nSignal);
% Compute a single metric for optimization purposes. Using simple mean
meanCorrScore = 0.5*(corrScoreX+corrScoreY);
corrScoreArray = [corrScoreX, corrScoreY];
end

%% Function used to compute objective for optimization
function [optScore, penaltyScore] = warpingObjective(optimWarp,nCtrlPts,inputSignals,nvArg)
% Control points are equally spaced in arc-length. 
% optimwarp is a column vector with first warped control point in the
% first nSignal indices, then 2nd control point in the next nSignal indices

% warpArray = reshape(optimWarp,length(inputSignals),nCtrlPts);
nSignal = length(inputSignals);
warpArray = reshape(optimWarp,[],nCtrlPts);
% Compute a warping penalty
penaltyScore = warpingPenalty(warpArray,nvArg.WarpingPenalty,nvArg);
penaltyScore = mean(penaltyScore);

% Perform warping - non-mex version
% [~, signalsX, signalsY] = warpArcLength(warpArray,inputSignals,...
%     nvArg.nResamplePoints);

% IMPORTANT: This is a compiled mex verison of warpArcLength. The mex
% function cannot be modified. If warpArcLength is updated later, you will
% also need to recompile the mex function
signalCellArray = cell(nSignal,1);
for i=1:nSignal
signalCellArray{i} = inputSignals(i).data;
end
[~, signalsX, signalsY] = warpArcLength_mex(warpArray,signalCellArray,nvArg.nResamplePoints);

% Compute correlation score
[corrScore, ~] = evalCorrScore(signalsX,signalsY);
% corrScore is a maximization goal. Turn into a minimization goal
optScore = 1-corrScore+penaltyScore;

end

%% Function used to warp arc-length
function [warpedSignals, signalsX, signalsY]...
    = warpArcLength(warpArray, inputSignals, nResamplePoints)
% Warp array: each row is warping points for an input signal, each column
% is warped point. Control points are interpolated  on [0,1] assuming
% equal spacing. 
nSignals = length(inputSignals);

% lmCtrlPts = linspace(0,1,2+nCtrlPts);
% lmCtrlPts = [0,warpArray(end,:),1];

% Initialize matrices
signalsX = zeros(nResamplePoints, nSignals);
signalsY = zeros(nResamplePoints, nSignals);
warpedSignals = cell(nSignals,1);

for iSignal = 1:nSignals
    % Assign responseSignal data array to matrix for brevity
    signal = inputSignals(iSignal).data;
    
    lmCtrlPts = [0,warpArray(iSignal+nSignals,:),1];
    
    % prepend 0 and append 1 to warp points for this signal to create valid
    % control points. 
    warpedCtrlPts = [0,warpArray(iSignal,:),1];
    
    % Construct warping function using SLM. This warps lmAlen to shiftAlen.
    % Use warping fuction to map computed arc-lengths onto the shifted
    % system. use built-in pchip function. This is a peicewise monotonic 
    % cubic spline. Signifincantly faster than SLM. 
    warpedNormAlen = pchip(lmCtrlPts,warpedCtrlPts,signal(:,4));
      
    % Now uniformly resample normalzied arc-length
    resamNormwarpedAlen = linspace(0,1, nResamplePoints)';
    resampX = interp1(warpedNormAlen, signal(:,1), resamNormwarpedAlen,'linear','extrap');
    resampY = interp1(warpedNormAlen, signal(:,2), resamNormwarpedAlen,'linear','extrap');
    % Assign to array for correlation calc
    signalsX(:,iSignal) = resampX;
    signalsY(:,iSignal) = resampY;
    
    % Assemble a cell array containing arrays of resampled signals. Similar
    % to 'normalizedSignal' in 'inputSignals' structure
    warpedSignals{iSignal} = [resamNormwarpedAlen,resampX,resampY];
end

end

%% Penalty function to prevent plateaus and extreme divergence in warping functions
function [penaltyScores] = warpingPenalty(warpArray,penaltyFactor,nvArg)
% Compute an array of penalty scores based on MSE between linear, unwarped
% arc-length and warped arc-length. Aim is to help prevent plateauing. 
[nSignals, nCtrlPts] = size(warpArray);
nSignals = nSignals/2;
% lmCtrlPts = [0, warpArray(end,:), 1];
penaltyScores = zeros(nSignals,1);
unwarpedAlen = linspace(0,1,nvArg.nResamplePoints);

for iSignal=1:nSignals
    penaltyScores(iSignal) = sum((unwarpedAlen - ...
        pchip([0,warpArray(iSignal+nSignals,:),1],...
        [0,warpArray(iSignal,:),1],unwarpedAlen)).^2);
end

penaltyScores = penaltyScores.*penaltyFactor;
end
    
