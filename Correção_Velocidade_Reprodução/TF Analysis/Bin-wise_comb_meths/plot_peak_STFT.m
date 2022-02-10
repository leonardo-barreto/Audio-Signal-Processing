Final_freq = fs/2;

freqs = 0 : Final_freq/length(SWGM_spectrogram(:,10)) : Final_freq - 1/length(SWGM_spectrogram(:,10));

figure;
hold on;
% plot(freqs, 10*log10(spectrograms_tensor(:,10,1)/max(spectrograms_tensor(:,10,1))), '--', 'linewidth', 2)
% plot(freqs, 10*log10(spectrograms_tensor(:,10,2)/max(spectrograms_tensor(:,10,2))), '--', 'linewidth', 2)
% plot(freqs, 10*log10(spectrograms_tensor(:,10,3)/max(spectrograms_tensor(:,10,3))), '--', 'linewidth', 2)
plot(freqs, 10*log10(spectrograms_tensor(:,10,1)), 'linewidth', 2)
plot(freqs, 10*log10(spectrograms_tensor(:,10,2)), 'linewidth', 2)
plot(freqs, 10*log10(spectrograms_tensor(:,10,3)), 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 2048', 'STFT 4096')

% axis([850 1150 -10 70])
xlabel ('Frequency [Hz]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_bin_wise_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);
