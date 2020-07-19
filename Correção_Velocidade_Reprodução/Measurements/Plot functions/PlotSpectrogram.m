function PlotSpectrogram(freqComponents,frameTimeInstants,powerMatrixDB)

    figure(1)
    surf(frameTimeInstants, freqComponents, powerMatrixDB)
    shading interp
    axis tight
    view(0, 90)
    set(gca, 'FontSize', 30,'yscale','log')
    xlabel('Tempo (s)','FontSize', 30)
    ylabel('Frequencia (Hz)','FontSize', 30)
    title('Espectrograma de potencia')

    hcol = colorbar;
    set(hcol, 'FontSize', 30)
    ylabel(hcol, 'Potencia (dB)')
    %ylim([0 3000])
end