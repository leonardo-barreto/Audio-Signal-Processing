function [spectrgMatrix,freqComponents,frameTimeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints)

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
        overlapSize = floor((overlapPerc/100)*windowSize); %This converts the overlap percentage to actual overlap size.
        [spectrgMatrix,f,frameTimeInstants,ps] = spectrogram (inputSignal,windowFunction,overlapSize,fftPoints,'power','onesided');
        freqComponents = transpose((samplingRate/2*pi).*f);
        powerMatrix = ps;
    else
        if DEBUG == 1
            fprintf(' Method used: stft function from Hristo Zhivomirov (MathWorks site).\n');
        end
        hopSize = floor(((100-overlapPerc)/100)*windowSize); %This converts the overlap percentage to hop size.
        [spectrgMatrix, freqComponents, frameTimeInstants] = stft(inputSignal,windowFunction,hopSize,fftPoints,samplingRate);
        powerMatrix = power(abs(spectrgMatrix),2)/fftPoints;
    end
    
end