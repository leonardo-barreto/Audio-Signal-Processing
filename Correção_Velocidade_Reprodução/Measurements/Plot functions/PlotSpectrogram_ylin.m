function PlotSpectrogram_ylin(freqComponents,frameTimes,plotLimits,powerSpectrgDB)

    figure;
    timeSize = floor(length(frameTimes/50));
    freqSize = floor(length(freqComponents/50));
    surf(frameTimes(1:timeSize), freqComponents(1:freqSize), powerSpectrgDB(1:freqSize,1:timeSize));

    %[t,f] = meshgrid(frameTimes, freqComponents);
    %surf(f,t,powerSpectrgDB);

    shading interp
    %axis tight
    view(0, 90)
    caxis(plotLimits)
    ylim([0 5000])
    
    set(gca, 'FontSize', 30)
    
    xlabel('Tempo (s)','FontSize', 30)
    ylabel('Frequencia (kHz)','FontSize', 30)
    %title('Espectrograma de potencia')

    hcol = colorbar;
    set(hcol, 'FontSize', 30)
    ylabel(hcol, 'Potencia (dB)')
    set(gca,'FontSize', 35)

end