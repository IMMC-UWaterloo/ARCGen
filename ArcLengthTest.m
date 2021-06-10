fclose all;
close all;
clear;
clc;

addpath('ThirdPartyFunctions')
stdFact = 1;
nResample = 100;

% % Dataset: Lessley (2004) Parabolas. 3  curves
% fNames = {'Lessley Parabolas/Dataset 0.csv',...
%     'Lessley Parabolas/Dataset 1.csv',...
%     'Lessley Parabolas/Dataset 2.csv',...
%     }
% index = [1,2];

% % Dataset: Lessley (2004) Parabolas with Outlier curve. 4 Curves
% fNames = {'Lessley Parabolas/Dataset 0.csv',...
%     'Lessley Parabolas/Dataset 1.csv',...
%     'Lessley Parabolas/Dataset 2.csv',...
%     'Lessley Parabolas/Dataset 3.csv',...
%     }
% index = [1,2];

% % Dataset: cRDCB Fixture Assessment 2 Force-Displacement
% fNames = {'cRDCB FA2 Force-Disp/FA2-01_PreFracture.csv',...
%     'cRDCB FA2 Force-Disp/FA2-02_PreFracture.csv',...
%     'cRDCB FA2 Force-Disp/FA2-03_PreFracture.csv',...
%     'cRDCB FA2 Force-Disp/FA2-04_PreFracture.csv',...
%     }
% index = [1,2];

% % Dataset: cRDCB Fixture Assessment 2 Traction-Separation
% fNames = {'cRDCB FA2 Trac-Sep/Processed_FA2-01_PreFractureSpline.csv',...
%     'cRDCB FA2 Trac-Sep/Processed_FA2-02_PreFractureSpline.csv',...
%     'cRDCB FA2 Trac-Sep/Processed_FA2-03_PreFractureSpline.csv',...
%     'cRDCB FA2 Trac-Sep/Processed_FA2-04_PreFractureSpline.csv',...
%     }
% index = [3,5];
% % index = [3,7];

% % Dataset: DCB 
% %   Force-Displacement: indices 2 & 3
% %   R-Curve: indices 4 & 5
% fNames = {...
%     'DCB Force-Disp/DCB-01_Toughness.csv',...
%     'DCB Force-Disp/DCB-02_Toughness.csv',...
%     'DCB Force-Disp/DCB-03_Toughness.csv',...
%     'DCB Force-Disp/DCB-05_Toughness.csv',...
%     'DCB Force-Disp/DCB-06_Toughness.csv',...
%     'DCB Force-Disp/DCB-07_Toughness.csv',...
%     'DCB Force-Disp/DCB-08_Toughness.csv',...
%     ...'DCB Force-Disp/DCB-09_Toughness.csv',...
%     }
% index = [2,3];
% % index = [4,5];

% % Watson 7333 Adhesive RDCB Force-Displacement Data
% fNames = {...
%     'Watson 7333 12R RDCB/12R1.csv',...
%     'Watson 7333 12R RDCB/12R3.csv',...
%     'Watson 7333 12R RDCB/12R4.csv',...
%     'Watson 7333 12R RDCB/12R6.csv',...
%     'Watson 7333 12R RDCB/12R7.csv',...
%     }
% index = [1,2];

% % Dataset: Kroell 1971 Thorax Impact Response, 16 MPH, 50 lb striker
% fNames =  {...
%     'Kroell 1971 Thorax Response/12FF.csv',...
%     'Kroell 1971 Thorax Response/13FM.csv',...
%     'Kroell 1971 Thorax Response/14FF.csv',...
%     'Kroell 1971 Thorax Response/15FM.csv',...
%     'Kroell 1971 Thorax Response/18FM.csv',...
%     'Kroell 1971 Thorax Response/19FM.csv',...
%     'Kroell 1971 Thorax Response/20FM.csv',...
%     'Kroell 1971 Thorax Response/22FM.csv',...
%     }
% index = [1,2]

% Mattucci Poastier Logitutinal Ligament Data - Quasi-static
fNames = {...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C080886_C2-C3.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C080886_C3-C5.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C080886_C6-C7.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C090033_C3-C4.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C090278_C3-C4.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C090278_C5-C6.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C090278_C7-T1.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C100923_C3-C4.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C100923_C5-C6.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/C100923_C7-T1.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/S090252_C7-T1.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/UB08L002_C5-C6.csv',...
    'Mattucci Ligament Data/Post Long Quasi - Pre-Failure/UB08L002_C7-T1.csv',...
    }
index = [1,2];
xPlotLim = [0,5];
yPlotLim = [0,750];

%% load data
for i=1:length(fNames)
    data(i).data = readmatrix(fNames{i});
    data(i).data = [data(i).data(:,index(1)),data(i).data(:,index(2))];
%     data(i).data = [...
%         smooth(data(i).data(:,index(1))-data(i).data(1,index(1)),10),...
%         smooth(data(i).data(:,index(2)),10)];
%     indexValid = and(data(i).data(:,1)<50, data(i).data(:,1)>0);
%     data(i).data = data(i).data(indexValid,:);
end

figure(); hold on;
for i=1:length(data)
    plot(data(i).data(:,1),data(i).data(:,2),'.-')
end
% Lobdell = readmatrix('Kroell 1971 Thorax Response/Lobdell 16mph Corridors.csv');
% xlim([0,4.5])
% ylim([0,1200])
% plot(Lobdell(:,1)+0.5,Lobdell(:,2)-150,'','LineWidth',3,'Color',[255, 213, 79]./255)
% plot(Lobdell(:,3)+0.5,Lobdell(:,4)-150,'o--','LineWidth',3,'Color',[255, 213, 79]./255)
% plot(Lobdell(:,5)+0.5,Lobdell(:,6)-150,'o--','LineWidth',3,'Color',[255, 213, 79]./255)

%% Compute arc-length
for i=1:length(data)
    temp = data(i).data;
    segments = sqrt( (temp(1:end-1,1)-temp(2:end,1)).^2 ...
        + (temp(1:end-1,2)-temp(2:end,2)).^2);
    alen = cumsum([0;segments]);
    data(i).data = [data(i).data,alen];
    data(i).maxAlen = alen(end);
    data(i).data = [data(i).data, alen./data(i).maxAlen];
    tempMax = max(temp,[],1);
    data(i).xMax = tempMax(1);
    data(i).yMax = tempMax(2);
end

fprintf('Average arc-length = %f +- %f\n',mean([data.maxAlen]),std([data.maxAlen]))
fprintf('Average x-max = %f +- %f\n',mean([data.xMax]),std([data.xMax]))
fprintf('Average y-max = %f +- %f\n',mean([data.yMax]),std([data.yMax]))


%% Resample curves based on normalized arc-length
for i=1:length(data)
	cfitx = fit(data(i).data(:,4),data(i).data(:,1),'linearinterp');
	cfity = fit(data(i).data(:,4),data(i).data(:,2),'linearinterp');
    
    normAlen = linspace(0,1,nResample)';
	data(i).normCurve = [normAlen,cfitx(normAlen),cfity(normAlen)];
end

figure();
subplot(2,2,[1,2]); hold on;
for i=1:length(data)
    plot(data(i).normCurve(:,2),data(i).normCurve(:,3),'.-')
end
xlabel('x-data')
ylabel('y-data')
title('Arc-length Discretized Normalized Curves')


%% Compute mean and std for each normalized arc-length
for i=1:nResample
    clear temp
    % collect each point
    for j=1:length(data)
        temp(j,:) = data(j).normCurve(i,2:3);
    end
    avgData(i,:) = mean(temp,1);
    stdevData(i,:) = std(temp,1);
end

subplot(2,2,[3]); hold on;
errorbar(data(1).normCurve(:,1),avgData(:,1),stdevData(:,1),'color',0.5.*[1,1,1])
for i=1:length(data)
    plot(data(i).normCurve(:,1),data(i).normCurve(:,2),'.-')
end
xlabel('Normalized Arc-length')
ylabel('x-data')
title('Average and St.Dev. of X-Data')

subplot(2,2,[4]); hold on;
errorbar(data(1).normCurve(:,1),avgData(:,2),stdevData(:,2),'color',0.5.*[1,1,1])
for i=1:length(data)
    plot(data(i).normCurve(:,1),data(i).normCurve(:,3),'.-')
end
xlabel('Normalized Arc-length')
ylabel('y-data')
title('Average and St.Dev. of Y-Data')

% %% Corridor development based on x,y data (not arc-length)
% figure(); hold on;
% plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',2.0)
% plot(avgData(:,1)+stdevData(:,1),avgData(:,2)+stdevData(:,2),'DisplayName','y+std')
% plot(avgData(:,1)-stdevData(:,1),avgData(:,2)+stdevData(:,2),'DisplayName','y-std')
% plot(avgData(:,1)+stdevData(:,1),avgData(:,2)-stdevData(:,2),'DisplayName','x+std')
% plot(avgData(:,1)-stdevData(:,1),avgData(:,2)-stdevData(:,2),'DisplayName','x-std')
% legend()
% 
% possibleCorridors = [avgData(:,1)+stdevData(:,1),avgData(:,2)+stdevData(:,2),...
%     avgData(:,1)-stdevData(:,1),avgData(:,2)+stdevData(:,2),...
%     avgData(:,1)+stdevData(:,1),avgData(:,2)-stdevData(:,2),...
%     avgData(:,1)-stdevData(:,1),avgData(:,2)-stdevData(:,2)];
% 
% % manual build corridors
% lowerCorr = [possibleCorridors(1:49,5:6);possibleCorridors(51:end,7:8)];
% upperCorr = [possibleCorridors(1:50,3:4);possibleCorridors(50:end,1:2)];
% 
% plot(lowerCorr(:,1),lowerCorr(:,2),'k--','LineWidth',1.5)
% plot(upperCorr(:,1),upperCorr(:,2),'k--','LineWidth',1.5)
% 
% figure(); hold on;
% plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',2.0)
% plot(lowerCorr(:,1),lowerCorr(:,2),'LineWidth',1.5,'DisplayName','Lower Corr.',...
%     'Color',[0.8,0.8,0.8])
% plot(upperCorr(:,1),upperCorr(:,2),'LineWidth',1.5,'DisplayName','Lower Corr.',...
%     'Color',[0.8,0.8,0.8])
% cmap = lines;
% for i=1:length(data)
%     plot(data(i).normCurve(:,2),data(i).normCurve(:,3),'.-',...
%         'DisplayName',['Data ' i],...
%         'Color', cmap(i,:))
% end
% 
% 
%% Corridor development using alighed ellipses
figure(); hold on;
plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',3.0)
% plot ellipses based on standard deviation
for i=1:nResample
    ellipse(stdevData(i,1).*stdFact, stdevData(i,2).*stdFact, 0,...
        avgData(i,1), avgData(i,2), 0.85.*[1,1,1])
end
cmap = lines;
for i=1:length(data)
    plot(data(i).data(:,1),data(i).data(:,2),'.-',...
        'DisplayName',['Data ' i],...
        'Color', cmap(i,:))
end
plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',2.0)
xlim(xPlotLim);
ylim(yPlotLim);

%% Corridor development using rectangles. Corners are +-st.dev.
figure(); hold on;
% plot ellipses based on standard deviation
for i=1:nResample
    rectangle('Position',...
        [avgData(i,1)-stdevData(i,1).*stdFact, avgData(i,2)-stdevData(i,2).*stdFact,...
        2.*stdevData(i,1).*stdFact, 2.*stdevData(i,2).*stdFact],...
        'EdgeColor',0.85*[1,1,1])
end
cmap = lines;
for i=1:length(data)
    plot(data(i).data(:,1),data(i).data(:,2),'.-',...
        'DisplayName',['Data ' i],...
        'Color', cmap(i,:))
end
plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',3.0)
xlim(xPlotLim);
ylim(yPlotLim);

% %% Corridor development using alighed ellipses, Kroell Specific
% cmap = cbrewer2('set1',5);
% figure(); hold on;
% plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',3.0)
% % plot ellipses based on standard deviation
% for i=1:nResample
%     ellipse(stdevData(i,1).*stdFact, stdevData(i,2).*stdFact, 0,...
%         avgData(i,1), avgData(i,2), 0.85.*[1,1,1])
% end
% for i=1:length(data)
%     plot(data(i).normCurve(:,2),data(i).normCurve(:,3),'-',...
%         'DisplayName',['Data ' i],...
%         'Color', brighten(cmap(2,:),0.7))
% end
% plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',2.0)
% 
% Lobdell = readmatrix('Kroell 1971 Thorax Response/Lobdell 16mph Corridors.csv');
% xlim([0,4.5])
% ylim([0,1200])
% plot(Lobdell(:,1)+0.5,Lobdell(:,2)-150,'','LineWidth',3,'Color',[255, 213, 79]./255)
% plot(Lobdell(:,3)+0.5,Lobdell(:,4)-150,'o--','LineWidth',3,'Color',[255, 213, 79]./255)
% plot(Lobdell(:,5)+0.5,Lobdell(:,6)-150,'o--','LineWidth',3,'Color',[255, 213, 79]./255)
% 
% %% Corridor development using rectangles. Corners are +-st.dev. Kroell Specific
% cmap = cbrewer2('set1',5);
% figure(); hold on;
% % plot ellipses based on standard deviation
% for i=1:nResample
%     rectangle('Position',...
%         [avgData(i,1)-stdevData(i,1).*stdFact, avgData(i,2)-stdevData(i,2).*stdFact,...
%         2.*stdevData(i,1).*stdFact, 2.*stdevData(i,2).*stdFact],...
%         'EdgeColor',0.85*[1,1,1])
% end
% for i=1:length(data)
%     plot(data(i).normCurve(:,2),data(i).normCurve(:,3),'.-',...
%         'DisplayName',['Data ' i],...
%         'Color', brighten(cmap(2,:),0.7))
% end
% plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',3.0)
% 
% Lobdell = readmatrix('Kroell 1971 Thorax Response/Lobdell 16mph Corridors.csv');
% xlim([0,4.5])
% ylim([0,1200])
% plot(Lobdell(:,1)+0.5,Lobdell(:,2)-150,'','LineWidth',3,'Color',[255, 213, 79]./255)
% plot(Lobdell(:,3)+0.5,Lobdell(:,4)-150,'o--','LineWidth',3,'Color',[255, 213, 79]./255)
% plot(Lobdell(:,5)+0.5,Lobdell(:,6)-150,'o--','LineWidth',3,'Color',[255, 213, 79]./255)