% Spectrograms combination
clc;
clear variables;
close all;

if isunix
    addpath ../Bin-wise_comb_meths
    addpath ./signals
    dirbar = '/';
else
    addpath ..\Bin-wise_comb_meths
    addpath .\signals
    dirbar = '\';
end

redLines = []; % plots red vertical dashed-lines for visual comparison of onsets
% sig_name = 'Ac_guit_note'; % for reading audio signal and naming the resulting figures
% sig_name = 'sin_1khz'; % for reading audio signal and naming the resulting figures
% sig_name = 'impulse'; % for reading audio signal and naming the resulting figures
sig_name = 'carnaval_percussao'; % for reading audio signal and naming the resulting figures

% sig_name = 'chirp_1_2000Hz_1s_sin_900Hz'; % for reading audio signal and naming the resulting figures
% alpha_ref = 1.9980; % alpha value at 0.5s to the signal chirp

plot_par = [0, .8, 0, 2000, 35]; % some parameters for plotting
plot_flags = [1, 1, 0];
fs = 44100/4;

figs_path = path_check(['.', dirbar, 'figs', dirbar, 'Lukin_Todd', dirbar, sig_name, dirbar]);

[data, fs_orig] = audioread([sig_name '.wav']);
if size(data,2) > 1
    x = mean(data.');
else
    x = data;
end
x = x(:);

% Resampling - - - - - - - -
if fs_orig ~= fs
    x = resample(x,fs,fs_orig);
end

% Standard setup
NFFT = 4096/4;
hop_DFT = 0.005; % time interval between frames in seconds

% % To plot 'sin_1khz'
% NFFT = 4096*2; % use this to plot 'sin_1khz'
% hop_DFT = 0.005; % time interval between frames in seconds

% % To plot 'impulse'
% NFFT = 4096/4;
% hop_DFT = 0.005; % time interval between frames in seconds

% Window lengths (in ascendent order!!!)
N_w = [1024 2048 4096]/4;
% N_w = [1024 4096]/4;
% alphas = [0 alpha_ref];
% alphas = [0 0 0 0];
samp_hop = ceil(hop_DFT*fs);
overlap_short = N_w(1) - samp_hop;
time_span = 0.05; % time span of the 2D window in seconds
Ncols_DFT = ceil(time_span/hop_DFT);
% Nrows_DFT = 10*N_w(end)/N_w(1);
Nrows_DFT = 10*NFFT/N_w(1);
% Nrows_DFT = 10*N_w(end)/N_w(1);
% size_w_S = [21 Ncols_DFT];
size_w_S = [Nrows_DFT  Ncols_DFT];

eta = 8;

%% Computing spectrograms - - - - - - -
[spectrograms_tensor, freq, t ] = spectrogram_tensor_prep(x, fs, N_w, NFFT, overlap_short);
% [spectrograms_tensor, freq, t ] = spectrogram_tensor_prep_FChT(x, fs, N_w, NFFT, alphas, overlap_short);

%% Combining spectrograms - - - - - -

% SWGM, just for comparison
% SWGM_spectrogram = SWGM_comb(spectrograms_tensor, 0.5);

combined_TFR_LT = spectrogram_comb_Lukin_Todd(spectrograms_tensor, size_w_S, eta);

computeGiniIndex(compress_dB_norm(combined_TFR_LT, 80))

% save('peak_LT')
% save('impulse_LT')

%% - - - - - - - - - 
% Showing images
% close all;

% plot_peak_SLS_LT
% plot_impulse_SLS_LT

% return

% N_w_str = '';
% for ind = 1:length(N_w)
%     N_w_str = [N_w_str sprintf('_%d', N_w(ind))];
%     if plot_flags(1)
%         plot_spec(sprintf('N = %d', N_w(ind)), compress_dB_norm(spectrograms_tensor(:,:,ind), 80), t, freq, plot_par, figs_path, redLines);
%     end
% end

if plot_flags(2)
plot_spec([sig_name '_LT'], compress_dB_norm(combined_TFR_LT, 80), t, freq, plot_par, figs_path, redLines);
%     plot_spec([sig_name '_SWGM'], compress_dB_norm(SWGM_spectrogram, 80), t, freq, plot_par, figs_path, redLines);
end                                                        