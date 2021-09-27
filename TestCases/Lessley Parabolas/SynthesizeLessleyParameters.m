fclose all;
close all;
clear;
clc;

coefs = [-240,380,0;
         -230,335,0;
         -245,310,0;
         -415,875,0]
specIds = {'Lessley Parabola 1',...
    'Lessley Parabola 2',...
    'Lessley Parabola 3',...
    'Lessley Parabola 4'}
     
figure(); hold on;

for i=1:size(coefs,1)
    paraRoots = roots(coefs(i,:));
    xx = linspace(0,paraRoots(2), 75)';
    yy = polyval(coefs(i,:),xx);
    plot(xx,yy,'.-');
    responseCurves(i).data = [xx,yy];
    responseCurves(i).specId = specIds{i}
end
