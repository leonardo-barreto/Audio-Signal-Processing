function [TFR_CQT, f, t] = TFAnalysis_CQT(x,fs)

    %   This function makes a time-frequency analysis of a signal, using the Constant-Q Transform.

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        hop = 128;
        N_w = 4096;
        NFFT = N_w;

    %% - - - - - - -  TFR computation - - - - - - -

        [spectrgMatrix, f, t] = spectrogram(x,hanning(N_w,'periodic'),N_w-hop,NFFT,fs);
        
        TFR_CQT = power(abs(spectrgMatrix),2);%/NFFT;
        %TFR_CQT = abs(spectrgMatrix);%/NFFT;
    
end