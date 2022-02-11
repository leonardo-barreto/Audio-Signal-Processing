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
        hop = 512;
        windowSize = 2048;
        fftPoints = 4096;

        % Compression
        plot_range = 80; %dB

    %% - - - - - - Input Reading - - - - - - 
        fprintf('\n\n------- READING INPUT ------\n\n')
        
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

    %% -------------------------------------- TFR stage -------------------------------------------
    fprintf('\n\n------- TIME-FREQUENCY ANALYSIS STARTED ------\n\n');
        [spectrgMatrix, freqComponents, timeInstants] = stft(x,hamming(windowSize,'periodic'),hop,fftPoints,fs);
        TFRepresentation = power(abs(spectrgMatrix),2)/fftPoints;
        
        fprintf(' Total frames: %i\n',length(timeInstants));
        fprintf(' Number of frequency bins: %i\n',length(freqComponents));
        fprintf('\nTime-Frequency Representation Finished.\n');

    %% ----------------- Plotting ---------------------
        
        PlotSpectrogram(freqComponents,timeInstants,10*log10(TFRepresentation));
    
    fprintf('\n\n------- TIME-FREQUENCY ANALYSIS FINISHED ------\n\n');
    
end