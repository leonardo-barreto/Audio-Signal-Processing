function y = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,DEBUG)

    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - STFT
    %
    %   2nd - Peak Detection
    %
    %   3rd - Frequency enhancements
    %
    %   4th - Tracking of sinusoidal components
    %


    % STFT and spectrogram stage

    if DEBUG == 1
        [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,1);
    else
        [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,0);
    end

    totalFrames = length(timeInstants); %Total number of signal frames
    totalFreqBins = length(freqComponents); %Total number of frequency components per frame

    % Calling Peak Detection

    if DEBUG == 1
        for frameCounter = 1:totalFrames
            detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(powerMatrix(:,frameCounter),totalFreqBins,samplingRate,windowSize,1);
        end
    else
        for frameCounter = 1:totalFrames
            detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(powerMatrix(:,frameCounter),totalFreqBins,samplingRate,windowSize,0);
        end
    end

    

end