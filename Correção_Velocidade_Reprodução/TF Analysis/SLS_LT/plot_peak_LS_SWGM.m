figure;
hold on;
plot(freq, 10*log10(spectrograms_tensor(:,100,1)/max(spectrograms_tensor(:,100,2))), 'linewidth', 2)
plot(freq, 10*log10(spectrograms_tensor(:,100,2)/max(spectrograms_tensor(:,100,2))), 'k', 'linewidth', 2)
plot(freq, 10*log10(SWGM_spectrogram(:,100)/max(spectrograms_tensor(:,100,2))), '--', 'linewidth', 2)
plot(freq, 10*log10(combined_TFR(:,100)/max(spectrograms_tensor(:,100,2))), '--', 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 4096', 'SWGM 0.5', 'LS')

axis([850 1150 -60 2])
xlabel ('Frequency [Hz]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_SWGM_LS'];
formatFig(gcf,figFileName,'en',figProp);