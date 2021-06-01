fclose all;
close all;
clear;
clc;

addpath('ThirdPartyFunctions')
% Dataset: Lessley (2004) Parabolas. 3  curves
% fNames = {'Lessley Parabolas/Dataset 0.csv',...
%     'Lessley Parabolas/Dataset 1.csv',...
%     'Lessley Parabolas/Dataset 2.csv'}

% Dataset: Lessley (2004) Parabolas with Outlier curve. 4 Curves
fNames = {'Lessley Parabolas/Dataset 0.csv',...
    'Lessley Parabolas/Dataset 1.csv',...
    'Lessley Parabolas/Dataset 2.csv',...
    'Lessley Parabolas/Dataset 3.csv'}

%% load data
for i=1:length(fNames)
    data(i).data = readmatrix(fNames{i});
end

figure(); hold on;
for i=1:length(data)
    plot(data(i).data(:,1),data(i).data(:,2),'.-')
end

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

%% Resample curves based on normalized arc-length
nResample = 100;
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
for i=1:length(data)
    plot(data(i).normCurve(:,1),data(i).normCurve(:,2),'.-')
end
errorbar(data(i).normCurve(:,1),avgData(:,1),stdevData(:,1),'color',0.0.*[1,1,1])
xlabel('Normalized Arc-length')
ylabel('x-data')
title('Average and St.Dev. of X-Data')

subplot(2,2,[4]); hold on;
for i=1:length(data)
    plot(data(i).normCurve(:,1),data(i).normCurve(:,3),'.-')
end
errorbar(data(i).normCurve(:,1),avgData(:,2),stdevData(:,2),'color',0.0.*[1,1,1])
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
plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',2.0)
% plot ellipses based on standard deviation
for i=1:nResample
    ellipse(stdevData(i,1), stdevData(i,2), 0,...
        avgData(i,1), avgData(i,2), 0.8.*[1,1,1])
end
cmap = lines;
for i=1:length(data)
    plot(data(i).normCurve(:,2),data(i).normCurve(:,3),'.-',...
        'DisplayName',['Data ' i],...
        'Color', cmap(i,:))
end

%% Corridor development using rectangles. Corners are +-st.dev.
figure(); hold on;
plot(avgData(:,1),avgData(:,2),'k','DisplayName','Char Avg','LineWidth',2.0)
% plot ellipses based on standard deviation
for i=1:nResample
    rectangle('Position',...
        [avgData(i,1)-stdevData(i,1), avgData(i,2)-stdevData(i,2),...
        2.*stdevData(i,1), 2.*stdevData(i,2)],...
        'EdgeColor',0.8*[1,1,1])
end
cmap = lines;
for i=1:length(data)
    plot(data(i).normCurve(:,2),data(i).normCurve(:,3),'.-',...
        'DisplayName',['Data ' i],...
        'Color', cmap(i,:))
end
