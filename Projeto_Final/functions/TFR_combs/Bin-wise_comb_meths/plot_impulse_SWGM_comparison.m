figure;
hold on;
plot(time, 10*log10(10^-10 + GM_spectrogram(1000, :)/max(GM_spectrogram(1000, :))), '--', 'linewidth', 2)
plot(time, 10*log10(10^-10 + SWGM_spectrogram(1000, :)/max(GM_spectrogram(1000, :))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + SWGM_spectrogram2(1000, :)/max(GM_spectrogram(1000, :))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + SWGM_spectrogram3(1000, :)/max(GM_spectrogram(1000, :))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + minimax_spectrogram(1000, :)/max(GM_spectrogram(1000, :))), '--', 'linewidth', 2)
hold off;

legend('GM', 'SWGM 0.1', 'SWGM 0.3', 'SWGM 0.5', 'MiniMax')

axis([.48 .52 -14 2])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_SWGM_comparison'];
formatFig(gcf,figFileName,'en',figProp);
