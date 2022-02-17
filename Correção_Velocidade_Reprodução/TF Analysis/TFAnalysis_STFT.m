function [freqComponents, timeInstants, TFRepresentation, spectrgMatrix] = TFAnalysis_STFT(x,fs)

    % This function computes a signal's STFT (by FFT).

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        %fs = 44100;
        hop = 128;
        windowSize = 2048;
        NFFT = windowSize;

        % Compression
        %plot_range = 80; %dB

    %% - - - - - - -  TFR computation - - - - - - -

        %[spectrgMatrix, freqComponents, timeInstants] = stft(x,hanning(windowSize,'periodic'),hop,NFFT,fs);
        [spectrgMatrix, freqComponents, timeInstants] = spectrogram(x,hanning(windowSize,'periodic'),windowSize-hop,NFFT,fs);
        
        TFRepresentation = power(abs(spectrgMatrix),2);%/NFFT;
        %TFRepresentation = abs(spectrgMatrix);%/NFFT;
        

    %% - - - - - - -  Plotting - - - - - - -  
        
        %PlotSpectrogram_ylin(freqComponents,timeInstants,[10-plot_range 10],10*log10(TFRepresentation));
    
end