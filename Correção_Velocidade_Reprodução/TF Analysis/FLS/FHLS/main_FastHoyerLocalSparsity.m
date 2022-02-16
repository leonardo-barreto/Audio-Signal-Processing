% Spectrograms combination
% clc;
% clear variables;
% close all;

if isunix
    addpath ../signals
    addpath ../util
    dirbar = '/';
else
    addpath ..\signals
    addpath ..\util
    dirbar = '\';
end

redLines = []; % plots red vertical dashed-lines for visual comparison of onsets
sig_name = 'opera_feminine_excerpt1'; % for reading audio signal and naming the resulting figures

resampling_factor = 1; % 1, 1/2 or 1/4

fs = 44100*resampling_factor;

[data, fs_orig] = audioread([sig_name '.wav']);

if size(data,2) > 1
    x = mean(data.');
else
    x = data;
end
x = x(:);

% Resampling
if fs_orig ~= fs
    x = resample(x,fs,fs_orig);
end

% % Adding noise
% x = x + 0.0001*randn(size(x));

sig_len = length(x)/fs;

%% Figures - - - - - - 
plot_par = [0, sig_len, 0, 2000, 35]; % some parameters for plotting
plot_flags = [0, 1]; % [input spectrograms, final spectrograms]

dir = ['.', dirbar, 'figs', dirbar, 'FastHoyerLocalSparsity', dirbar, sig_name, dirbar];

figs_path = path_check(dir);

%% System Setup

% Standard setup
NFFT = 4096*resampling_factor;
hop_DFT = 0.005; % Final time interval between frames in seconds (interp after all computations)
% hop_DFT = 0.001; % Final time interval between frames in seconds (interp after all computations)

% % To plot 'sin_1khz'
% NFFT = 4096*2; % use this to plot 'sin_1khz'

% % To plot 'impulse'
% NFFT = 4096/4;
 
% % Window lengths
% N_w = [1, 2, 3, 4]*1024*resampling_factor;
% % N_w = [1024 4096]/4;
% overlap_short = N_w(1)*0.75; % overlap to be used initially in the spectrograms
% % overlap_percentage = 0.75; % overlap to be used initially in the spectrograms
% hop_DFT = 0.005; % Final time interval between frames in seconds (interp after all computations)
% % hop_DFT = 0.001; % Final time interval between frames in seconds (interp after all computations)
% % final_time_scale = 0:hop_DFT:length(x)/fs;
% 
% % 2D Analysis window
% time_span = 0.05; % time span of the 2D window in seconds
% freq_span = 1722; % frequency span of the 2D window in Hz
% size_W_t_f = [freq_span , time_span]; % In terms of Hz and s

% Window lengths (in ascendent order!!!)
N_w = [1, 2, 3, 4]*1024*resampling_factor;
samp_hop = ceil(hop_DFT*fs);
overlap_short = N_w(1) - samp_hop;

% 2D Analysis window
time_span = 0.05; % time span of the 2D window in seconds
Ncols_DFT = ceil(time_span/hop_DFT);
Nrows_DFT = 3*N_w(end)/N_w(1);
size_W_m_k = [Nrows_DFT Ncols_DFT];

% size_w_E = [Nrows_DFT/3   2*Ncols_DFT]; %%%%%%%%%%

% plot_window_E(size_w_E, 'window_E.png', 40)
% plot_window_S(size_w_S, 'window_S.png', 40)
% 
% return

eta = 20;

%% Computing spectrograms - - - - - - -
[spectrograms_tensor, f, t] = spectrogram_tensor_prep(x, fs, N_w, NFFT, overlap_short);

combined_TFR_HLS = spectrogram_comb_FastHoyerLocalSparsity(spectrograms_tensor, size_W_m_k, eta);


%% Showing images - - - - - - - - - 
% close all;

N_w_str = '';
for ind = 1:length(N_w)
% for ind = [1, 3]
    N_w_str = [N_w_str sprintf('_%d', N_w(ind))];
    if plot_flags(1)
        plot_spec(sprintf('N = %d', N_w(ind)), compress_dB_norm(spectrograms_tensor(:,:,ind), 80), ...
                                            t, f, plot_par, figs_path, redLines);
    end
%     fprintf('Energy: %d\n', sum(sum(spectrograms_tensor(:,:,ind))))
end

if plot_flags(2)
    plot_spec([sig_name '_FHLS'], compress_dB_norm(combined_TFR_HLS, 80), ...
        t, f, plot_par, figs_path, redLines);
end                                                        


% if plot_flags(2)
%     plot_spec([sig_name '_SLS'], compress_dB_norm(combined_TFR_SLS, 80), ...
%         t, f, plot_par, figs_path, redLines);
% end                                                        