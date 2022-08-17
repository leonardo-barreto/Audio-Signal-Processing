% Spectrograms combination
clc;
clear variables;
close all;

if isunix
    dirbar = '/';
else
    dirbar = '\';
end

% signal_name = 'sin_1khz'; % for reading audio signal and naming the resulting figures
% signal_name = 'impulse'; % for reading audio signal and naming the resulting figures
% signal_name = 'impulse_sin'; % for reading audio signal and naming the resulting figures
% signal_name = 'Ac_guit_note'; % for reading audio signal and naming the resulting figures
%signal_name = 'Paulistana1'; % for reading audio signal and naming the resulting figures

% some parameters for plotting [t_i, t_f, f_min, f_max, font_size]
% plot_par = [0, 0.2, 0, 4000, 25]; 
plot_par = [0, 1, 0, 5000, 35]; 
plot_flags = [1, 1]; % 1- Input TFRs, 2- Combined TFR

reduction_factor = 1;
fs = 48000/reduction_factor;

% Windows length
N_w = [1024 2048 4096]/reduction_factor;
hop = 128/reduction_factor;
N_FFT = 4096/reduction_factor;

overlap_short = N_w(1) - hop;

[data, fs_orig] = audioread([signal_name '.wav']);
if size(data,2) > 1
    x = mean(data.');
else
    x = data;
end
x = x(:);
x = x(1:fs);

% % Resampling - - - - - - - -
% x = resample(x,fs,fs_orig);

% Plot signal and spectrogram - - - - - -
figs_path = path_check(['.', dirbar, 'figs', dirbar, signal_name, dirbar]);

[spec, freq, time] = spectrogram(x, 2048, 2048 - hop, N_FFT, fs);

plot_spec(sprintf('%s', signal_name), compress_dB_norm(abs(spec),80), time, freq, plot_par, figs_path);
return
% Plot signal and spectrogram - - - - - -



%% Computing TFRs

[spectrograms_tensor, freq, time] = spectrogram_tensor_prep(x, fs, N_w, overlap_short, N_FFT);

% % Combining representations - - - - - -
% SWGM_spectrogram = SWGM_comb(spectrograms_tensor, 0.1);
% SWGM_spectrogram2 = SWGM_comb(spectrograms_tensor, 0.3);
 SWGM_spectrogram3 = SWGM_comb(spectrograms_tensor, 0.5);
% GM_spectrogram = SWGM_comb(spectrograms_tensor, 0);
% minimax_spectrogram = minimax_comb(spectrograms_tensor);
% NM_spectrogram = NM_comb(spectrograms_tensor);
% RM_spectrogram = RM_comb(spectrograms_tensor);

%% - - - - - - - - - - - - - - - - - -
% Showing images
% close all;

figs_path = path_check(['.', dirbar, 'figs', dirbar, signal_name, dirbar]);

% plot_peak_NM
% plot_peak_RM
% plot_peak_GM
% plot_peak_minimax
 plot_peak_SWGM
% % plot_peak_SWGM_comparison % MUST REMOVE ENERGY COMPENSATION
% plot_peak_bw_methods

% plot_impulse
% plot_impulse_SWGM_comparison % MUST REMOVE ENERGY COMPENSATION

% N_w_str = '';
% for ind = 1:length(N_w)
%     N_w_str = [N_w_str sprintf('_%d', N_w(ind))];
%     if plot_flags(1)
%         plot_spec(sprintf('%s_%d', signal_name, N_w(ind)), spectrograms_tensor(:,:,ind).^.25, time, freq, plot_par, figs_path);
%     end
% end
% 
% % N_w_str = '';
% % for ind = 1:length(N_w)
% %     N_w_str = [N_w_str sprintf('_%d', N_w(ind))];
% %     if plot_flags(1)
% %         plot_spec(sprintf('N = %d', N_w(ind)), spectrograms_tensor(:,:,ind).^.25, time, freq, plot_par, figs_path);
% %     end
% % end
% 
% fig_name = sprintf('SWGM Product(Beta=%1.1f)%s', beta, N_w_str);
% fig_name(fig_name == '.') = '_';
% if plot_flags(2)
%     plot_spec(fig_name, SWGM_spectrogram, time, freq, plot_par, figs_path);
% end

% disp('All the figures are now available in directory "figs".')

% close all;