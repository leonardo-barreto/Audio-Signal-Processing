function PlotPVC(variationCurve,timeInstants)
    
    %{
    deviationPercentageMatrix = (deviationMatrix-1).*100; %This for plotting purposes only.  
    figure;
    surf(1:totalFrames,1:totalTracks,deviationPercentageMatrix);
    shading interp
    axis tight
    view(0, 90)
    xlabel('Frames')
    ylabel('Tracks (Indexes)')
    title('Deviation curve of each track around its own mean frequency');
    hcol = colorbar;
    ylabel(hcol, 'Deviation from mean frequency (%)');
    %}

    figure;
    plot(timeInstants,variationCurve,'LineWidth',3);
    X = sprintf('Curva de desvio relativo de pitch');
    title(X,'FontSize', 30);
    xlabel('Tempo (s)','FontSize', 30);
    ylabel('Frequencia relativa','FontSize', 30);
    set(gca,'FontSize', 30)
    ylim([min(variationCurve) max(variationCurve)])
end