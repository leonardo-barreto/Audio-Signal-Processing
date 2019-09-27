function y = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,DEBUG)

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

    %if DEBUG == 1
        %[spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,1);
    %else
        [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc);
    %end

    powerMatrixDB = 10*log(powerMatrix);

    totalFreqBins = length(freqComponents);
    totalFrames = length(timeInstants);

    %Building a signal frame as a structure with its fields
    signalFrame = {};
    signalFrame.totalFrames = totalFrames; %Total number of signal frames
    signalFrame.currentFrame = 1;
    signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
    signalFrame.freqComponents = freqComponents; %frequency components vector


    %Random frame chosen for plots (temporary)
    DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);

    % Calling Peak Detection

    detectedPeaksMatrix = zeros(totalFreqBins,totalFrames);

    if DEBUG == 1
        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,samplingRate,1,DEBUG_FRAME);
        end
    else
        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,samplingRate,0,0);
        end
    end

    %TEMPORARY
    y = detectedPeaksMatrix;

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

            %if (initialPeaks(freqCounter))
             %   stem(freqComponents(freqCounter),initialPeaks(freqCounter));
            %elseif (detectedPeaks(freqCounter))
            %    stem(freqComponents(freqCounter),detectedPeaks(freqCounter),'filled');
            %else
            %    plot(freqComponents(freqCounter),powerSpectrumDB(freqCounter));
            %end
        %end

        X = sprintf('Peak Detection of frame %i of %i',DEBUG_FRAME,totalFrames);
        title(X);
        xlabel('Frequency (kHz)');
        ylabel('Power (dB)');
        legend ('Power Spectrum','Detected Peaks');

        hold off;
    end

    

end