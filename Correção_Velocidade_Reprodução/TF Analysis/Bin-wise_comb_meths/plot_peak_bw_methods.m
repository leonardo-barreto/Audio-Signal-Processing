Final_freq = fs/2;

freqs = 0 : Final_freq/length(SWGM_spectrogram(:,10)) : Final_freq - 1/length(SWGM_spectrogram(:,10));

norm_GM_spectrogram = GM_spectrogram/sqrt(sum(sum(GM_spectrogram.^2)));
norm_minimax_spectrogram = minimax_spectrogram/sqrt(sum(sum(minimax_spectrogram.^2)));
norm_SWGM_spectrogram3 = SWGM_spectrogram2/sqrt(sum(sum(SWGM_spectrogram2.^2)));
norm_NM_spectrogram = NM_spectrogram/sqrt(sum(sum(NM_spectrogram.^2)));
norm_RM_spectrogram = RM_spectrogram/sqrt(sum(sum(RM_spectrogram.^2)));

figure;
hold on;
plot(freqs, 10*log10(norm_NM_spectrogram(:,10)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(freqs, 10*log10(norm_RM_spectrogram(:,10)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(freqs, 10*log10(norm_GM_spectrogram(:,10)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(freqs, 10*log10(norm_minimax_spectrogram(:,10)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(freqs, 10*log10(norm_SWGM_spectrogram3(:,10)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
hold off;

legend('NM', 'RM', 'GM', 'MiniMax', 'SWGM 0.5')

axis([850 1150 -60 2])
xlabel ('Frequency [Hz]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_bin_wise_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);

axis([980 1020 -5 1])

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_bin_wise_comb_comparison_zoom'];
formatFig(gcf,figFileName,'en',figProp);
