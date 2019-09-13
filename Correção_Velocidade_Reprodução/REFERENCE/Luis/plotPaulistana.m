clc; close all; clear;
load paulistana.mat
load paulistana2.mat
timeAxis = centralSamples/fs;
timeAxis2 = centralSamples2/fs;

figure;

subplot(2,1,1);
plot(timeAxis,smoothCurve,'linewidth',1.5);
ylabel('PVC','interpreter','latex');
xlabel('Time (s)','interpreter','latex');
ylim([0.993 1.007]);
xlim([1 4])

subplot(2,1,2); 
plot(timeAxis2,movmean(smoothCurve2,20),'linewidth',1.5);
ylabel('PVC','interpreter','latex');
xlabel('Time (s)','interpreter','latex');
ylim([0.993 1.007]);
xlim([1 4])

% figProp = struct('size',12,'font','Times','lineWidth',1,'figDim',[1 1 600 350]);
% formatFig(gcf,'paulistana','en',figProp);