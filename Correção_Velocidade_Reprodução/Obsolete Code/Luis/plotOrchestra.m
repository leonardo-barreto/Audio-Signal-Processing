clc; close all; clear;
load orchestra.mat
load orchestra2.mat
timeAxis = centralSamples/fs;
timeAxis2 = centralSamples2/fs;
t = 1:nFrames;
A = 0.022;
phase = 1.8;

f = 1;
q = 178;

true = 1 + A*cos(2*pi*f*t/q + phase);

figure;

subplot(2,1,1);
plot(timeAxis,smoothCurve,'linewidth',1.5);
hold on;
plot(timeAxis,true,'k:','linewidth',1.5);
ylabel('PVC','interpreter','latex');
xlabel('Time (s)','interpreter','latex');
ylim([0.97 1.03]);
xlim([min(timeAxis) max(timeAxis)])

subplot(2,1,2); 
plot(timeAxis2,smoothCurve2,'linewidth',1.5);
ylabel('PVC','interpreter','latex');
xlabel('Time (s)','interpreter','latex');
ylim([0.97 1.03]);
xlim([min(timeAxis2) max(timeAxis2)])

figProp = struct('size',12,'font','Times','lineWidth',1,'figDim',[1 1 600 350]);
formatFig(gcf,'orchestra','en',figProp);