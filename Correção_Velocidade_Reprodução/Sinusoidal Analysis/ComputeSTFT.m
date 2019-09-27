function [spectrgMatrix,freqComponents,frameTimeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc)

    % This function computes a signal's STFT (by FFT). It includes windowing and spectrogram generation.
    
    fftPoints = windowSize; %DFT is calculated through all points of the window.

    overlapSize = (overlapPerc/100)*windowSize; %This converts the overlap percentage to actual overlap size.

    switch windowType
        case 'hann' 
            windowFunction = hann(windowSize);
        case 'hamming'
            windowFunction = hamming(windowSize);
        otherwise
            error('Invalid window type. Valid options are ''hann'' or ''hamming''.');
    end

    [s,f,t,ps] = spectrogram (inputSignal,windowFunction,overlapSize,fftPoints,samplingRate,'power','onesided');

    spectrgMatrix = s;
    freqComponents = f;
    frameTimeInstants = t;
    powerMatrix = ps;

end