function [f, t, combined_TFR_FLS] = TFAnalysis_FLS(x,fs)
    
    %   This function makes a time-frequency analysis of a signal, using the FLS method.

    %% - - - - - - Parameters configuration - - - - - - 

        % Input Signal params
        resampling_factor = 1; % 1, 1/2 or 1/4

        % Standard setup
        NFFT = 4096*resampling_factor;
        hop_DFT = 0.0029; % Final time interval between frames in seconds (interp after all computations)
        % hop_DFT = 0.001; % Final time interval between frames in seconds (interp after all computations)

        % Window lengths (in ascendent order!!!)
        N_w = [1, 2, 3, 4]*1024*resampling_factor;
        %samp_hop = ceil(hop_DFT*fs);
        samp_hop = 128;
        overlap_short = N_w(1) - samp_hop;

        % 2D Analysis window
        time_span = 0.03; % time span of the 2D window in seconds
        Ncols_DFT = ceil(time_span/hop_DFT);
        Nrows_DFT = 3*N_w(end)/N_w(1);
        size_W_m_k = [Nrows_DFT Ncols_DFT];

        gamma = 20;

    %% - - - - - - Input Reading - - - - - - 

        

    %% - - - - - - -  TFR computation - - - - - - - 
 
        % - - - - - - - Computing FLS-combined spectrograms - - - - - - -
        [spectrograms_tensor, f, t] = spectrogram_tensor_prep(x, fs, N_w, NFFT, overlap_short);

        combined_TFR_FLS = spectrogram_comb_FastHoyerLocalSparsity(spectrograms_tensor, size_W_m_k, gamma);
                 
        %combined_TFR_FLS = compress_dB_norm(combined_TFR_FLS, plot_range);
        %combined_TFR_FLS = 10*log10(combined_TFR_FLS);    

end
