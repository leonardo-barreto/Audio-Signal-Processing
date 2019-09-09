function PlotTonalness(CurTonalness,curPeaks,fftLen,fs,CurMagSpecData)

th=0.8;
xdb = 20*log10(CurMagSpecData);
fAxis = (0:fftLen/2)*fs/fftLen;
fig = figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);
yyaxis left
plot(fAxis,CurTonalness,'m','lineWidth',2);
hold on;
plot(fAxis(curPeaks),CurTonalness(curPeaks), 'r*','MarkerSize',10);
hold on
plot([fAxis(1) fAxis(end)],[th th],'k--');

xlabel('Frequency (Hz)','interpreter','latex');
ylabel('$\mathcal{T}(k,b)$','interpreter','latex');
ylim([-0.05 1.05]);
xlim([0 2000]);
set(gca,'YTick',[0 0.5 1]);

yyaxis right;
plot(fAxis,xdb,'k');
ylabel('$|X(k,b)|$ (dB)','interpreter','latex');
ylim([-33 65]);
set(gca,'YTick',[-20 0 20 40 60]);
figProp = struct('size',12,'font','Times','lineWidth',1,'figDim',[1 1 600*0.85 380*0.85]);
formatFig(gcf,'figTonalness','pt',figProp);

end