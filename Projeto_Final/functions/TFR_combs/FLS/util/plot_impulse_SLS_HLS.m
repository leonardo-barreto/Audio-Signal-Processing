load('impulse_SLS')

figure;
hold on;
plot(t, 10*log10(10^-10 + spectrograms_tensor(200, :, 1)/max(spectrograms_tensor(200, :, 1))), 'k', 'linewidth', 2)
plot(t, 10*log10(10^-10 + spectrograms_tensor(200, :, 2)/max(spectrograms_tensor(200, :, 1))), 'linewidth', 2)
plot(t, 10*log10(10^-10 + combined_TFR(200, :)/max(spectrograms_tensor(200, :, 1))), '--', 'linewidth', 2)
load('impulse_HLS')
plot(final_time_scale, 10*log10(10^-10 + combined_TFR_HLS(200, :)/max(max(spectrograms{1,1}))), '--', 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 4096', 'SLS', 'HLS')

axis([.45 .55 -30 2])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_SLS_HLS'];
formatFig(gcf,figFileName,'en',figProp);
