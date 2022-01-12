function PlotPeakDetection(sinAnalysisParameters,signalFrame)

    if (isempty(signalFrame.peakMatrix))
        fprintf('\nNo peaks in this frame. Well, nothing to plot.\n\n');
        return
    end

    figure
    hold on;

    plot(signalFrame.freqComponents./1000,signalFrame.powerSpectrumDB,'G','LineWidth',2);
    plot(signalFrame.freqComponents./1000,signalFrame.powerSpectrumThreshold,'B','LineWidth',2);


    plot(signalFrame.peakMatrix(2,:)/1000,signalFrame.peakMatrix(1,:),'r.','markersize',30);

    X = sprintf('Quadro %i de %i',signalFrame.currentFrame,sinAnalysisParameters.totalFrames);
        title(X,'FontSize', 30);
        xlabel('Frequencia (kHz)','FontSize', 30);
        ylabel('Potencia (dB)','FontSize', 30);
        legend ('Espectro','Limiar SSE','Picos detectados');
        set(gca,'FontSize', 30);
        xlim([0 signalFrame.freqComponents(end)/1000])
        hold off;

end