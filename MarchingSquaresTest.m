close all;
clear;
clc;

addpath('Functions')

% Array format [x0, y0, rx, ry] one row per ellipse

ellipses = [0,0,1,2;
            0.7,1.7,2,1;
            -0.5,1.5,1,1;
            -2,-1,1.5,2;
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

%% start sorting algorithm
% find the minimum y-value
[min,ind] = min(vertices(:,2));

% initalize first point
envelop = vertices(ind,1:2);
% add second point from same array index
envelop = [envelop; vertices(ind,3:4)];
% indexUsed = zeros(size(vertices,1),1);
% indexUsed(ind) = 1;
indicesUsed = ind;

for i = 2:2*size(vertices,1)

% find the repeat of the 2nd point
col1 = find(all(vertices(:,1:2)==repmat(envelop(end,:),size(vertices,1),1),2))
if size(col1,1)>1
    indexExtra = find(any(col1'==indicesUsed,2));
    indexExtra = find(col1~=indicesUsed(indexExtra));
    col1 = col1(indexExtra);
end
    
col2 = find(all(vertices(:,3:4)==repmat(envelop(end,:),size(vertices,1),1),2))
if size(col2,1)>1
    indexExtra = find(any(col2'==indicesUsed,2));
    indexExtra = find(col2~=indicesUsed(indexExtra));
    col2 = col2(indexExtra);
end

if isempty(col1)
    col1 = indicesUsed(end)
end
if isempty(col2)
    col2 = indicesUsed(end)
end

if any(repmat(col1,size(indicesUsed,1),1) == indicesUsed)
    envelop = [envelop;vertices(col2,1:2)];
    indicesUsed = [indicesUsed;col2];
elseif any(repmat(col2,size(indicesUsed,1),1)==indicesUsed)
    envelop = [envelop;vertices(col1,3:4)];
    indicesUsed = [indicesUsed;col1];
end

end


figure(); hold on
comet(envelop(:,1),envelop(:,2))
plot(envelop(:,1),envelop(:,2),'x-','MarkerSize',12)
scatter([vertices(:,1);vertices(:,3)],[vertices(:,2);vertices(:,4)])

    

