function [spectrgMatrix,powerMatrix,freqComponents,frameTimeInstants] = TFAnalysis_STFT(inputSignal,samplingRate,windowType,windowSize,hopSize,fftPoints)

    % This function computes a signal's STFT (by FFT).

    METHOD = 1;
    DEBUG = 0;

    switch windowType
        case 'hann' 
            windowFunction = hann(windowSize,'periodic');
        case 'hamming'
            windowFunction = hamming(windowSize,'periodic');
        case 'rectangular'
            windowFunction = rectwin(windowSize);
        case 'triangular'
            windowFunction = triang(windowSize);
        otherwise
            error('Invalid window type. Valid options are ''hann'', ''hamming'', ''rectangular'' or ''triangular''.');
    end

    if METHOD == 0
        if DEBUG == 1
            fprintf(' Method used: MATLAB spectrogram.\n');
        end
        [spectrgMatrix,f,frameTimeInstants,ps] = spectrogram (inputSignal,windowFunction,hopSize,fftPoints,'power','onesided');
        freqComponents = transpose((samplingRate/2*pi).*f);
        powerMatrix = ps;
    else
        if DEBUG == 1
            fprintf(' Method used: stft function from Hristo Zhivomirov (MathWorks site).\n');
        end
        [spectrgMatrix, freqComponents, frameTimeInstants] = stft(inputSignal,windowFunction,hopSize,fftPoints,samplingRate);
        powerMatrix = power(abs(spectrgMatrix),2)/fftPoints;
    end
    
end