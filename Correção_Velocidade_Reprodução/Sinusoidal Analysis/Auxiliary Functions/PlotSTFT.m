function PlotSTFT(S,f,t)

    figure(1)
    surf(t, f, S)
    shading interp
    axis tight
    view(0, 90)
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
    xlabel('Time, s')
    ylabel('Frequency, Hz')
    title('Amplitude spectrogram of the signal')

    hcol = colorbar;
    set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
    ylabel(hcol, 'Magnitude, dB')
    
end