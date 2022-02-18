function [freqComponents, timeInstants, TFRepresentation, spectrgMatrix] = TFAnalysis_STFT(x,fs)

    % This function computes a signal's STFT (by FFT).

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        hop = 128;
        %N_w = 2570;
        N_w = 4096;
        NFFT = N_w;

    %% - - - - - - -  TFR computation - - - - - - -

        %[spectrgMatrix, freqComponents, timeInstants] = stft(x,hanning(N_w,'periodic'),hop,NFFT,fs);
        [spectrgMatrix, freqComponents, timeInstants] = spectrogram(x,hanning(N_w,'periodic'),N_w-hop,NFFT,fs);
        
        TFRepresentation = power(abs(spectrgMatrix),2);%/NFFT;
        %TFRepresentation = abs(spectrgMatrix);%/NFFT;
    
end