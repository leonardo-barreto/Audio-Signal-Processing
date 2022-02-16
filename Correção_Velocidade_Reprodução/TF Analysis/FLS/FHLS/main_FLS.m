% Spectrograms combination
% clc;
% clear variables;
% close all;

% - - - - - - Setting Paths - - - - - - 

    if isunix
        figsPath = path_check('./audio/figs/');
        dirbar = '/';
    else
        figsPath = path_check('.\audio\figs\');
        dirbar = '\';
    end

    if isunix
        addpath ./FLS/util
        dirbar = '/';
    else
        addpath .\FLS\util
        dirbar = '\';
    end

%% - - - - - - - Plotting parameters - - - - - - - 
    %redLines = []; % plots red vertical dashed-lines for visual comparison of onsets

%% - - - - - - - Figure parameters - - - - - - 
    %sig_len = length(x)/fs;
    %plot_par = [0, sig_len, 0, 2000, 35]; % some parameters for plotting
    %plot_flags = [0, 0]; % [input spectrograms, final spectrograms]

    %dir = ['.', dirbar, 'audio', dirbar, 'figs', dirbar, 'FLS', dirbar, signal_name, dirbar];

    %figs_path = path_check(dir);

%% - - - - - - - Computing FLS-combined spectrograms - - - - - - -
    [spectrograms_tensor, f, t] = spectrogram_tensor_prep(x, fs, N_w, NFFT, overlap_short);

    combined_TFR_HLS = spectrogram_comb_FastHoyerLocalSparsity(spectrograms_tensor, size_W_m_k, eta);


%% - - - - - - - Plotting - - - - - - - - - 

    
%{
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
            plot_spec([signal_name '_FHLS'], compress_dB_norm(combined_TFR_HLS, 80), ...
                t, f, plot_par, figs_path, redLines);
        end         
%}
                                                                                                   