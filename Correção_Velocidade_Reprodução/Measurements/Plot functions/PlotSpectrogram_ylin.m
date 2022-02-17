function PlotSpectrogram_ylin(freqs,times,plotLimits,powerSpectrgDB)

    figure;

    %% SURF
        
        %{
        timeSize = floor(length(times/50));
        freqSize = floor(length(freqs/50));
        surf(times(1:timeSize), freqs(1:freqSize), powerSpectrgDB(1:freqSize,1:timeSize), 'EdgeColor', 'none');

        %axis tight
        view(0, 90)
        %}
        
    %% IMAGESC

        imagesc(times, freqs, 10*log10(powerSpectrgDB));
        set(gca,'YDir','normal');
        colormap(1-gray);
    
    %% SHARED AND LABELS

        caxis(plotLimits)
        ylim([0 10000])

        xlabel('Tempo (s)','FontSize', 30)
        ylabel('Frequencia (kHz)','FontSize', 30)
        
        %hcol = colorbar;
        %set(hcol, 'FontSize', 30)
        %ylabel(hcol, 'Potencia (dB)')
        set(gca,'FontSize', 15)

end