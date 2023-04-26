function [combined_TFR_MRFCI, f, t] = TFAnalysis_MRFCI(x,fs)
    
    %   This function makes a time-frequency analysis of a signal, using the MRFCI method.

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -
        hop = 128;
        block_size = 2; % Analysis block size in seconds

        % Synthesis TFR - - - - - - - - - - - - - - - - - - -
        N_alphas = 7; % Number of alphas (I, in the paper)
        Nf = [1024 2048 3072 4096]; % Window sizes
        asym_flags = [0 0 0 1]; % Indicate which Nf will be used to compute Fan-chirp instances
        FChT_flags = [0 0 1 1]; % Indicate which Nf will be used to compute Fan-chirp instances

        % Structure tensor - - - - - - - - -
        C_limits = [0 .25 .75 1 ; 1 2 2 3]; % Defines the points of transition for the lambda weights
        range = 100; % range in dB for the ST

        % G parameters
        sigma_t = Nf(end)/(4*fs); % in s
        sigma_f = 100;% in Hz

        Nf_structure_tensor = Nf(2); % To compute the structure tensor

    %% - - - - - - -  TFR computation - - - - - - -

        % Computing combination procedure in blocks
        [combined_TFR_MRFCI, STFT_tensor, t, f] = MRFCI(x, fs, Nf, block_size,...
        FChT_flags, asym_flags, Nf_structure_tensor, hop, N_alphas, C_limits, range, sigma_t, sigma_f);
            
        %Transposing because MRFCI function outputs 1 x N vector
        f = f';
    

end
