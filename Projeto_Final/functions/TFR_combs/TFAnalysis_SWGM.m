function [combined_TFR_SWGM, f, t] = TFAnalysis_SWGM(x,fs)
    
    %   This function makes a time-frequency analysis of a signal, using the SWGM method. 

    %% - - - - - - Parameters configuration - - - - - -

        N_w = [1024 2048 4096 8192];
        hop = 128;
        N_FFT = 8192;
        overlap_short = N_w(1) - hop;

        beta = 0.5; % Amount of weight applied

    %% - - - - - - -  TFR computation - - - - - - -

        [spectrograms_tensor, f, t] = spectrogram_tensor_prep_binwise(x, fs, N_w, overlap_short, N_FFT);
        combined_TFR_SWGM = SWGM_comb(spectrograms_tensor, beta);

end
