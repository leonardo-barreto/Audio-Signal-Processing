Final_freq = fs/2;

freqs = 0 : Final_freq/length(SWGM_spectrogram(:,10)) : Final_freq - 1/length(SWGM_spectrogram(:,10));

figure;
hold on;
plot(freqs, 10*log10(GM_spectrogram(:,10)/max(GM_spectrogram(:,10))), '--', 'linewidth', 2)
plot(freqs, 10*log10(SWGM_spectrogram(:,10)/max(GM_spectrogram(:,10))), 'linewidth', 2)
plot(freqs, 10*log10(SWGM_spectrogram2(:,10)/max(GM_spectrogram(:,10))), 'linewidth', 2)
plot(freqs, 10*log10(SWGM_spectrogram3(:,10)/max(GM_spectrogram(:,10))), 'linewidth', 2)
plot(freqs, 10*log10(minimax_spectrogram(:,10)/max(GM_spectrogram(:,10))), '--', 'linewidth', 2)
hold off;

legend('GM', 'SWGM 0.1', 'SWGM 0.3', 'SWGM 0.5', 'MiniMax')

axis([850 1150 -60 2])
xlabel ('Frequency [Hz]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_SWGM_comparison'];
formatFig(gcf,figFileName,'en',figProp);
