function [freqComponents,timeInstants,TFRepresentation,spectrgMatrix] = TFAnalysis_STFT(signal_name)

    % This function computes a signal's STFT (by FFT).

    if isunix
        addpath ./audio
        dirbar = '/';
    else
        addpath .\audio
        dirbar = '\';
    end
    
    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        fs = 44100;
        hop = 128;
        windowSize = 4096;
        NFFT = 4096;

        % Compression
        plot_range = 80; %dB

    %% - - - - - - Input Reading - - - - - - 
        
        [data, fs_orig] = audioread([signal_name]);

        if size(data,2) > 1
            x = mean(data.');
        else
            x = data;
        end
        x = x(:);

        if fs ~= fs_orig
            x = resample(x, fs, fs_orig);
        end

    %% - - - - - - -  TFR computation - - - - - - -

        %[spectrgMatrix, freqComponents, timeInstants] = stft(x,hanning(windowSize,'periodic'),hop,NFFT,fs);
        [spectrgMatrix, freqComponents, timeInstants] = spectrogram(x,hanning(windowSize,'periodic'),windowSize-hop,NFFT,fs);
        
        TFRepresentation = power(abs(spectrgMatrix),2);%/NFFT;
        %TFRepresentation = abs(spectrgMatrix);%/NFFT;
        

    %% - - - - - - -  Plotting - - - - - - -  
        
        %PlotSpectrogram_ylin(freqComponents,timeInstants,[10-plot_range 10],10*log10(TFRepresentation));
    
end