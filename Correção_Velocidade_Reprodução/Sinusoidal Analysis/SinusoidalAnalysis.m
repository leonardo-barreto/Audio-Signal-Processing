function y = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,DEBUG)

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

    fprintf('\nSinusoidal Analysis started.\n Sampling Rate(Hz): %i\n Window: %s (size %i, overlap %i%%) \n FFT Points: %i\n', samplingRate,windowType,windowSize,overlapPerc,fftPoints);

    [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);
    
    powerMatrixDB = 10*log(powerMatrix);

    totalFreqBins = length(freqComponents);
    totalFrames = length(timeInstants);

    %Building a signal frame as a structure with its fields
    signalFrame = {};
    signalFrame.totalFrames = totalFrames; %Total number of signal frames
    signalFrame.currentFrame = 1; %Current frame
    signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
    signalFrame.freqComponents = freqComponents; %frequency components vector


    % Calling Peak Detection

    fprintf('\nPeak Detection Started.\n');

    detectedPeaksMatrix = zeros(totalFreqBins,totalFrames);

    if DEBUG == 1
        %Random frame chosen for DEBUG (temporary)
        DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);

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
        detectedPeakPositions = find(detectedPeaksMatrix(:,DEBUG_FRAME));

        figure
        plot(freqComponents,powerSpectrumDB,'B');
        hold on;
        
        for freqCounter = detectedPeakPositions
            plot(freqComponents(freqCounter),detectedPeaks(freqCounter),'r*');
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

    amplitudeMatrixDB = 10*log(abs(spectrgMatrix));
    

end