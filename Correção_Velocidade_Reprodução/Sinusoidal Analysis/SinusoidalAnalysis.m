function y = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,thresholdMethod,DEBUG)

    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - STFT
    %
    %   2nd - Frequency Enhancement
    %
    %   3rd - Peak Detection
    %
    %   4th - Tracking of sinusoidal components
    %


    % STFT and spectrogram stage

    fprintf('\n\nSinusoidal Analysis started.\n Sampling Rate: %i\n Window: %s of size %i and %i%% overlap\n FFT Points: %i\n Threshold Method: %s\n', samplingRate,windowType,windowSize,overlapPerc,fftPoints,thresholdMethod);

    [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);
    
    powerMatrixDB = 10*log(powerMatrix);

    totalFreqBins = length(freqComponents);
    totalFrames = length(timeInstants);

    %Building a signal frame as a structure with its fields
    signalFrame = {};
    signalFrame.totalFrames = totalFrames; %Total number of signal frames
    signalFrame.currentFrame = 1;
    signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
    signalFrame.freqComponents = freqComponents; %frequency components vector


    %Random frame chosen for DEBUG (temporary)
    if DEBUG == 1
        DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);
    end

    % Calling Peak Detection

    detectedPeaksMatrix = zeros(totalFreqBins,totalFrames);

    if (DEBUG == 1)
        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            if DEBUG_FRAME == signalFrame.currentFrame
                detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,1);
            else
                detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,0);
            end
        end
    else
        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,0);
        end
    end

    if DEBUG == 1
        powerSpectrumDB = powerMatrixDB(:,DEBUG_FRAME);
        detectedPeaks = detectedPeaksMatrix(:,DEBUG_FRAME);

        figure
        plot(freqComponents,powerSpectrumDB,'B');
        hold on;
        for freqCounter = 1:totalFreqBins
            if (detectedPeaks(freqCounter))
                plot(freqComponents(freqCounter),detectedPeaks(freqCounter),'r*');
            end
        end

        X = sprintf('Peak Detection of frame %i of %i',DEBUG_FRAME,totalFrames);
        title(X);
        xlabel('Frequency (kHz)');
        ylabel('Power (dB)');
        legend ('Power Spectrum','Detected Peaks');

        hold off;
    end

    %TEMPORARY
    y = detectedPeaksMatrix;
    

end