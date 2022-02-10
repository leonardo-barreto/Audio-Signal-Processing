% Spectrograms combination
clc;
clear variables;
close all;

if isunix
    dirbar = '/';
else
    dirbar = '\';
end

range = [40 70 100];
axis_par = [0.3, 0.7, 800, 1200];
plot_flags = [1, 1]; % 1- Input TFRs, 2- Combined TFR
resample_flag = 1;
fs = 44100/4;
% Windows length
N_w = [1024 2048 4096]/4;
hop = 64/4;
overlap_short = N_w(1) - hop;
sig_name = 'impulse_sin'; % for reading audio signal and naming the resulting figures
figs_path = path_check(['.', dirbar, 'figs', dirbar, sig_name, dirbar]);

[impulse, fs_orig] = audioread('impulse.wav');
[sin_1khz, ~] = audioread('sin_1khz.wav');

data = 512*impulse + sin_1khz + randn(size(sin_1khz))*10^(0);

if size(data,2) > 1
    x = mean(data.');
else
    x = data;
end
x = x(:);

% Resampling
if resample_flag
    x = resample(x,fs,fs_orig);
else
    fs = fs_orig;
end

%%
[spectrograms_tensor, freq, time] = spectrogram_tensor_prep(x, fs, N_w, overlap_short);

% Combining representations - - - - - -
% SWGM_spectrogram = SWGM_comb(spectrograms_tensor, 0.1);
% SWGM_spectrogram2 = SWGM_comb(spectrograms_tensor, 0.3);
SWGM_spectrogram3 = SWGM_comb(spectrograms_tensor, 0.5);
GM_spectrogram = SWGM_comb(spectrograms_tensor, 0);
minimax_spectrogram = minimax_comb(spectrograms_tensor);
NM_spectrogram = NM_comb(spectrograms_tensor);
RM_spectrogram = RM_comb(spectrograms_tensor);

%% Plot all bin-wise combinations side-by-side  - - - - - - - - - - -

figure;
% NM
subplot(5,3,1)
imagesc(time, freq, compress_dB_norm(NM_spectrogram, range(1)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
% ylabel ('Frequency [Hz]');
title('NM')

subplot(5,3,2)
imagesc(time, freq, compress_dB_norm(NM_spectrogram, range(2)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('NM')

subplot(5,3,3)
imagesc(time, freq, compress_dB_norm(NM_spectrogram, range(3)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('NM')

% RM
subplot(5,3,4)
imagesc(time, freq, compress_dB_norm(RM_spectrogram, range(1)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
% ylabel ('Frequency [Hz]');
title('RM')

subplot(5,3,5)
imagesc(time, freq, compress_dB_norm(RM_spectrogram, range(2)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('RM')

subplot(5,3,6)
imagesc(time, freq, compress_dB_norm(RM_spectrogram, range(3)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('RM')

% GM
subplot(5,3,7)
imagesc(time, freq, compress_dB_norm(GM_spectrogram, range(1)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
% ylabel ('Frequency [Hz]');
title('GM')

subplot(5,3,8)
imagesc(time, freq, compress_dB_norm(GM_spectrogram, range(2)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('GM')

subplot(5,3,9)
imagesc(time, freq, compress_dB_norm(GM_spectrogram, range(3)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('GM')

% MM
subplot(5,3,10)
imagesc(time, freq, compress_dB_norm(minimax_spectrogram, range(1)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
% ylabel ('Frequency [Hz]');
title('MM')

subplot(5,3,11)
imagesc(time, freq, compress_dB_norm(minimax_spectrogram, range(2)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('MM')

subplot(5,3,12)
imagesc(time, freq, compress_dB_norm(minimax_spectrogram, range(3)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
title('MM')

% SWGM
subplot(5,3,13)
imagesc(time, freq, compress_dB_norm(SWGM_spectrogram3, range(1)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
% xlabel ('Time [s]');
% ylabel ('Frequency [Hz]');
title('SWGM')

subplot(5,3,14)
imagesc(time, freq, compress_dB_norm(SWGM_spectrogram3, range(2)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
% xlabel ('Time [s]');
title('SWGM')

subplot(5,3,15)
imagesc(time, freq, compress_dB_norm(SWGM_spectrogram3, range(3)))
colormap(abs(1-gray))
axis(axis_par);
set(gca,'YDir','normal');
% xlabel ('Time [s]');
title('SWGM')

figProp = struct('size', 20, 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 800 2000]);
figFileName = [figs_path sig_name '_bin_wise_comb_comparison_TFR'];
formatFig(gcf, figFileName, 'en', figProp);

% save(figFileName, 'tit', 'spectro', 't', 'freq', 'plot_par');