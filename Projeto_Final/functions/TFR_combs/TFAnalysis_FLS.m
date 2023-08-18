function [combined_TFR_FLS, f, t] = TFAnalysis_FLS(x,fs)
    
    %   This function makes a time-frequency analysis of a signal, using the FLS method. 

    %% - - - - - - Parameters configuration - - - - - -

        NFFT = 8192;
        hop_DFT = 0.0029; % Final time interval between frames in seconds (interp after all computations)
        N_w = [1, 2, 3, 4]*1024;
        hop = 128;
        overlap_short = N_w(1) - hop;
        % 2D Analysis window
        time_span = 0.03; % time span of the 2D window in seconds
        Ncols_DFT = ceil(time_span/hop_DFT);
        Nrows_DFT = 3*N_w(end)/N_w(1);
        size_W_m_k = [Nrows_DFT Ncols_DFT];
        % Contrast factor for weights
        gamma = 20;

    %% - - - - - - -  TFR computation - - - - - - -

        [spectrograms_tensor, f, t] = spectrogram_tensor_prep(x, fs, N_w, NFFT, overlap_short);
        combined_TFR_FLS = spectrogram_comb_FastHoyerLocalSparsity(spectrograms_tensor, size_W_m_k, gamma);

end
