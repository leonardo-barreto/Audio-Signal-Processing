function PlotSpectrogram_tiled(freqs,times,plotLimits,powerSpectrgDB)

    figure;
    tiledlayout()

    %% SURF
        
        %{
        timeSize = floor(length(times/50));
        freqSize = floor(length(freqs/50));
        surf(times(1:timeSize), freqs(1:freqSize), powerSpectrgDB(1:freqSize,1:timeSize), 'EdgeColor', 'none');

        %axis tight
        view(0, 90)
        %}
        
    %% IMAGESC

        imagesc(times, freqs, powerSpectrgDB);
        set(gca,'YDir','normal');
        colormap(1-gray);
    
    %% SHARED AND LABELS

        caxis(plotLimits)
        ylim([0 20000])

        xlabel('Tempo (s)','FontSize', 15)
        ylabel('Frequencia (kHz)','FontSize', 15)
        
        hcol = colorbar;
        set(hcol, 'FontSize', 15)
        ylabel(hcol, 'Potencia (dB)')
        set(gca,'FontSize', 15)

end