function PlotSpectrogram(freqComponents,frameTimeInstants,powerMatrixDB)

    figure;
    timeSize = floor(length(frameTimeInstants/50));
    freqSize = floor(length(freqComponents/50));
    surf(frameTimeInstants(1:timeSize), freqComponents(1:freqSize), powerMatrixDB(1:freqSize,1:timeSize));
    shading interp
    axis tight
    view(0, 90)

    %Log scale
    %set(gca, 'FontSize', 30, 'yscale', 'log')
    %Linear scale
    set(gca, 'FontSize', 30)

    ylim([0 5000])

    xlabel('Tempo (s)','FontSize', 30)
    ylabel('Frequencia (kHz)','FontSize', 30)
    %title('Espectrograma de potencia')

    hcol = colorbar;
    set(hcol, 'FontSize', 30)
    ylabel(hcol, 'Potencia (dB)')
    set(gca,'FontSize', 35)
    %ylim([0 3000])
end