function y = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,DEBUG)

    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - STFT
    %
    %   2nd - Frequency Enhancement
    %
    %   3rd - Peak Detection
    %
    %   4th - Construction of sinusoidal tracks
    %


    % -------------------------------------- STFT and spectrogram stage -------------------------------------------

        fprintf('\nSinusoidal Analysis started.\n Sampling Rate(Hz): %i\n Window: %s (size %i, overlap %i%%) \n FFT Points: %i\n', samplingRate,windowType,windowSize,overlapPerc,fftPoints);

        [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);
        
        powerMatrixDB = 10*log(powerMatrix); 

        totalFreqBins = length(freqComponents);
        totalFrames = length(timeInstants);

        % Building a signal frame as a peak detection entity
        signalFrame = {};
        signalFrame.totalFrames = totalFrames; %Total number of signal frames
        signalFrame.currentFrame = 1; %Current frame
        signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
        signalFrame.freqComponents = freqComponents; %frequency components vector
        signalFrame.samplingRate = samplingRate;
        signalFrame.fftPoints = fftPoints;

    % ---------------------------------------------- Peak Detection ------------------------------------------------

        fprintf('\nPeak Detection Started.\n');

        detectedPeaksMatrix = {};
        detectedFrequenciesMatrix = {};
        peakMatrix = {};

        if DEBUG == 1
            %Random frame chosen for DEBUG (temporary)
            %DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);

            DEBUG_FRAME = 392;

            for frameCounter = 1:totalFrames
                signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
                signalFrame.currentFrame = frameCounter;
                if DEBUG_FRAME == signalFrame.currentFrame
                    [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,1);
                else
                    [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,0);
                end
            end

        else

            for frameCounter = 1:totalFrames
                signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
                signalFrame.currentFrame = frameCounter;
                [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate);
            end

        end

        for index = 1:totalFrames
            peakMatrix{index} = [detectedPeaksMatrix{index};detectedFrequenciesMatrix{index}];
        end

    y = peakMatrix;

end