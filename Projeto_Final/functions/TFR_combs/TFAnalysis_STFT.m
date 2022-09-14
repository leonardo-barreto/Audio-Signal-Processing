function [TFR_STFT, f, t] = TFAnalysis_STFT(x,fs)

    %   This function makes a time-frequency analysis of a signal, using the traditional STFT method.

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        hop = 128;
        N_w = 4096;
        NFFT = N_w;

    %% - - - - - - -  TFR computation - - - - - - -

        [spectrgMatrix, f, t] = spectrogram(x,hanning(N_w,'periodic'),N_w-hop,NFFT,fs);
        
        TFR_STFT = power(abs(spectrgMatrix),2);%/NFFT;
        %TFR_STFT = abs(spectrgMatrix);%/NFFT;
    
end