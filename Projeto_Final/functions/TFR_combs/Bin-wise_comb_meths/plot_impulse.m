% Plot NM

figure;
hold on;
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 1)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 2)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 3)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + NM_spectrogram(1000, :)/max(spectrograms_tensor(1000, :, 1))), '--k', 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 2048', 'STFT 4096', 'Num. Mean')

axis([.45 .55 -30 5])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_NM_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);

% Plot RM

figure;
hold on;
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 1)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 2)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 3)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + RM_spectrogram(1000,:)/max(spectrograms_tensor(1000, :, 1))), '--k', 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 2048', 'STFT 4096', 'Recip. Mean')

axis([.45 .55 -30 5])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_RM_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);

% Plot GM 

figure;
hold on;
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 1)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 2)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 3)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + GM_spectrogram(1000, :)/max(spectrograms_tensor(1000, :, 1))), '--k', 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 2048', 'STFT 4096', 'Geom. Mean')

axis([.45 .55 -30 5])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_GM_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);

% Plot minimax

figure;
hold on;
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 1)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 2)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 3)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + minimax_spectrogram(1000, :)/max(spectrograms_tensor(1000, :, 1))), '--k', 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 2048', 'STFT 4096', 'Minimax')

axis([.45 .55 -30 5])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_minimax_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);

% Plot SWGM

figure;
hold on;
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 1)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 2)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + spectrograms_tensor(1000, :, 3)/max(spectrograms_tensor(1000, :, 1))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + SWGM_spectrogram3(1000, :)/max(spectrograms_tensor(1000, :, 1))), '--k', 'linewidth', 2)
hold off;

legend('STFT 1024', 'STFT 2048', 'STFT 4096', 'SWGM 0.5')

axis([.45 .55 -30 5])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_SWGM_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);

% Plot all bin-wise combination methods

norm_GM_spectrogram = GM_spectrogram/sqrt(sum(sum(GM_spectrogram.^2)));
norm_minimax_spectrogram = minimax_spectrogram/sqrt(sum(sum(minimax_spectrogram.^2)));
norm_SWGM_spectrogram3 = SWGM_spectrogram2/sqrt(sum(sum(SWGM_spectrogram2.^2)));
norm_NM_spectrogram = NM_spectrogram/sqrt(sum(sum(NM_spectrogram.^2)));
norm_RM_spectrogram = RM_spectrogram/sqrt(sum(sum(RM_spectrogram.^2)));

figure;
hold on;
plot(time, 10*log10(10^-10 + norm_NM_spectrogram(1000, :)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + norm_RM_spectrogram(1000, :)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + norm_GM_spectrogram(1000, :)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + norm_minimax_spectrogram(1000, :)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
plot(time, 10*log10(10^-10 + norm_SWGM_spectrogram3(1000, :)/max(max(norm_SWGM_spectrogram3))), 'linewidth', 2)
hold off;

legend('NM', 'RM', 'GM', 'MiniMax', 'SWGM 0.5')

axis([.45 .55 -30 5])
xlabel ('Time [s]');
ylabel ('Magnitude [dB]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path sig_name '_bin_wise_comb_comparison'];
formatFig(gcf,figFileName,'en',figProp);
