close all;
clear;
clc;

addpath('Functions')

% Array format [x0, y0, rx, ry] one row per ellipse

ellipses = [0,0,1,2;
            0.7,1.7,2,1;
            -0.5,1.5,1,1;
            -2,-1,1.5,2;
            1.65, 0.75, 0.5, 1.5;
            ];

xlims = [-4,4];
ylims = [-4,4];

nPts = 50;

[xx,yy] = meshgrid(linspace(xlims(1),xlims(2),nPts),...
    linspace(ylims(1),ylims(2),nPts));

zz = zeros(size(xx));

for iPt = 1:nPts
    for jPt = 1:nPts
        for kEllip = 1:size(ellipses,1)
            zz(iPt,jPt) = max(zz(iPt,jPt),...
                ((xx(iPt,jPt) - ellipses(kEllip,1)).^2 ./ ellipses(kEllip,3).^2 + ...
                (yy(iPt,jPt) - ellipses(kEllip,2)).^2 ./ ellipses(kEllip,4).^2).^-1);
        end
    end
end
        
figure()
contourf(xx,yy,zz,[1,1])
colorbar()
hold on
axis equal

cmap = cbrewer2('set2',2);

%% Start marching squares
figure()
hold on
for iEllip = 1:size(ellipses,1)
    ellipse(ellipses(iEllip,3),ellipses(iEllip,4),0,ellipses(iEllip,1),ellipses(iEllip,2))
end

colormap(cmap);
pScat = scatter(xx(:),yy(:),20,zz(:)>=1,'Filled');
axis equal

vertices = [];
for iPt = 1:(nPts-1)  % Rows (y-axis)
    for jPt = 1:(nPts-1)   % Columns (x-axis)
        % Cell values
        %  1 -- 2 
        %  |    |
        %  |    |
        %  8 -- 4
        %
        % REMEMBER, DUMBASS: 
        % array(i,j) = array(rows, columns,) = array(y,x)

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
                vertices = [vertices;
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]];
                
                line([interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)) , xx(iPt,jPt) ],...
                    [yy(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)) ],...
                    'color',[0,0,0],'LineWidth',2);         
            case 3
                % West-North
                vertices = [vertices;
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                
                line([xx(iPt+1,jPt), interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    [interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)), yy(iPt+1,jPt)],...
                    'color',[0,0,0],'LineWidth',2);
            case 4
                % North-South
                vertices = [vertices;
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt) ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                
                line([interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)), ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1))],...
                    [yy(iPt,jPt), yy(iPt+1,jPt)],...
                    'color',[0,0,0],'LineWidth',2);
            case 5
                % North-East
                vertices = [vertices;...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];
                
                line([interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),xx(iPt+1,jPt+1)],...
                    [yy(iPt+1,jPt+1), interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2);           
            case 6  % Ambiguous
                % South-West
                vertices = [vertices;...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]];
                
                line([interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)) , xx(iPt,jPt)],...
                    [yy(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)) ],...
                    'color',[0,0,0],'LineWidth',2); 
                
                % North-East
                vertices = [vertices;...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];
                
                line([interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),xx(iPt+1,jPt+1)],...
                    [yy(iPt+1,jPt+1), interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2);      
            case 7
                % West-East
                vertices = [vertices;...
                    [xx(iPt,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
                
                line([xx(iPt,jPt), xx(iPt,jPt+1)],...
                    [interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2);  
            case 8
                % South - East
                vertices = [vertices;...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
                
                line([interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),xx(iPt,jPt+1)],...
                    [yy(iPt,jPt+1), interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2);    
            case 9
                % South - East
                vertices = [vertices;...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
                
                line([interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),xx(iPt,jPt+1)],...
                    [yy(iPt,jPt+1), interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2);
            case 10
                % West-East
                vertices = [vertices;...
                    [xx(iPt,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
                
                line([xx(iPt,jPt), xx(iPt,jPt+1)],...
                    [interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2);  
            case 11 % Ambiguous
                % West-North
                vertices = [vertices;
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                
                line([xx(iPt+1,jPt), interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    [interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)), yy(iPt+1,jPt)],...
                    'color',[0,0,0],'LineWidth',2);
                
                % South-East
                vertices = [vertices;...
                    [interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),yy(iPt,jPt+1),...
                    xx(iPt,jPt+1),interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))]];
                
                line([interpVal(xx(iPt,jPt+1),zz(iPt,jPt+1),xx(iPt,jPt),zz(iPt,jPt)),xx(iPt,jPt+1)],...
                    [yy(iPt,jPt+1), interpVal(yy(iPt,jPt+1),zz(iPt,jPt+1),yy(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2);
            case 12
                % North-East
                vertices = [vertices;...
                    [interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt+1),...
                    xx(iPt+1,jPt+1),interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))]];
                
                line([interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),xx(iPt+1,jPt+1)],...
                    [yy(iPt+1,jPt+1), interpVal(yy(iPt+1,jPt+1),zz(iPt+1,jPt+1), yy(iPt,jPt+1), zz(iPt,jPt+1))],...
                    'color',[0,0,0],'LineWidth',2); 
            case 13
                % North-South
                vertices = [vertices;
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt) ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                
                line([interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)), ...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt),xx(iPt+1,jPt+1), zz(iPt+1,jPt+1))],...
                    [yy(iPt,jPt), yy(iPt+1,jPt)],...
                    'color',[0,0,0],'LineWidth',2);
            case 14
                % West-North
                vertices = [vertices;
                    [xx(iPt+1,jPt),interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)),...
                    interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1)),yy(iPt+1,jPt)]];
                
                line([xx(iPt+1,jPt), interpVal(xx(iPt+1,jPt),zz(iPt+1,jPt), xx(iPt+1,jPt+1),zz(iPt+1,jPt+1))],...
                    [interpVal(yy(iPt,jPt),zz(iPt,jPt), yy(iPt+1,jPt),zz(iPt+1,jPt)), yy(iPt+1,jPt)],...
                    'color',[0,0,0],'LineWidth',2);
            case 15
                % South-West
                vertices = [vertices;...
                    [interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)),yy(iPt,jPt),...
                    xx(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt))]];
                
                line([interpVal(xx(iPt,jPt),zz(iPt,jPt),xx(iPt,jPt+1),zz(iPt,jPt+1)) , xx(iPt,jPt) ],...
                    [yy(iPt,jPt), interpVal(yy(iPt,jPt),zz(iPt,jPt),yy(iPt+1,jPt),zz(iPt+1,jPt)) ],...
                    'color',[0,0,0],'LineWidth',2);   
            case 16
                % No vertices 
        end
    end
end

scatter([vertices(:,1);vertices(:,3)],[vertices(:,2);vertices(:,4)])

% comet(vertices(:,1),vertices(:,2))

figure(); hold on
for iEllip = 1:size(ellipses,1)
    ellipse(ellipses(iEllip,3),ellipses(iEllip,4),0,ellipses(iEllip,1),ellipses(iEllip,2))
end

%% start sorting algorithm
% find the minimum y-value
[min,lastIndex] = min(vertices(:,2));

% initalize first point
envelope = vertices(lastIndex,1:2);
% add second point from same array index
envelope = [envelope; vertices(lastIndex,3:4)];
indexUsed = zeros(size(vertices,1),1);
indexUsed(lastIndex) = 1;

exitFlag = 0;   % Set exit flag fudge
% Go though all vertices looking for the next connecting face
for iVerts = 2:size(vertices,1)
    % For an enclosed polygon, all lines share vertices
    % Find the all repeated vertices
    foundVert12 = find(all(ismembertol(vertices(:,1:2), envelope(end,:)), 2));
    foundVert34 = find(all(ismembertol(vertices(:,3:4), envelope(end,:)), 2));
    
    % there will only ever be two points, distributed between foundVert12
    % and foundVert34. Select the vertex which is NOT the same as the last
    % vertex. This indices a new line segments. 
    if ~isempty(foundVert12)
        for iInd = 1:length(foundVert12)
            if foundVert12(iInd) ~= lastIndex
                envelope = [envelope; vertices(foundVert12(iInd),3:4)];
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
                envelope = [envelope; vertices(foundVert34(iInd),1:2)];
                lastIndex = foundVert34(iInd);
                indexUsed(lastIndex) = 1;
                break;
            end
        end
    end
    
end
% plot(envelop(:,1),envelop(:,2),'x-','MarkerSize',12)
comet(envelope(:,1),envelope(:,2))
line1 = [-3, -3; -2, -2];
line2 = [2, 2; 3, 3];
line(line1(:,1),line1(:,2))
line(line2(:,1),line2(:,2))

%Find intercepts to divide line using Poly
[xInt,yInt,iIntLow] = polyxpoly(envelope(:,1),envelope(:,2),line1(:,1),line1(:,2))
plot(xInt,yInt,'*k')
plot(envelope(iIntLow,1),envelope(iIntLow,2),'sg')
iIntLow = iIntLow(1);

[xInt,yInt,iIntUpp] = polyxpoly(envelope(:,1),envelope(:,2),line2(:,1),line2(:,2))
plot(xInt,yInt,'*k')
plot(envelope(iIntUpp,1),envelope(iIntUpp,2),'sg')
iIntUpp = iIntUpp(1);

% To find inner or outer corridors, first determine if polygon is clockwise
% or counter-clockwise. Then, based on which index is large, separate out
% inner and outer corridor based on which intercept index is larger. 
if ispolycw(envelope(:,1),envelope(:,2))
    if iIntLow > iIntUpp
        outerCorr = [envelope(iIntLow:end,:);envelope(1:iIntUpp,:)];
        innerCorr = envelope(iIntUpp:iIntLow,:);
    else
        outerCorr = envelope(iIntLow:iIntUpp,:);
        innerCorr = [envelope(iIntUpp:end,:);envelope(1:iIntLow,:)];
    end
else
    if iIntLow > iIntUpp
        innerCorr = [envelope(iIntLow:end,:);envelope(1:iIntUpp,:)];
        outerCorr = envelope(iIntUpp:iIntLow,:);
    else
        innerCorr = envelope(iIntLow:iIntUpp,:);
        outerCorr = [envelope(iIntUpp:end,:);envelope(1:iIntLow,:)];
    end
end

plot(innerCorr(:,1),innerCorr(:,2),'o-r')
plot(outerCorr(:,1),outerCorr(:,2),'o-m')


