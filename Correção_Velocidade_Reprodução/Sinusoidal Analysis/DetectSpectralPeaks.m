function [detectedFinalPeaks,detectedPeakFrequencies] = DetectSpectralPeaks(inputFrame,samplingRate,DEBUG)

    % This function aims to detect spectral peaks in a signal's frame, given its spectrum.
    %
    % Firstly, a threshold estimation method is applied.
    %
    % Secondly, the peaks are detected based on the threshold.
    %
    %

    % Gathering frame data
    
        powerSpectrumDB = inputFrame.powerSpectrumDB;
        freqComponents = inputFrame.freqComponents;
        totalFreqBins = inputFrame.totalFreqBins;
        totalFrames = inputFrame.totalFrames;
        currentFrame = inputFrame.currentFrame;
        samplingRate = inputFrame.samplingRate;
        fftPoints = inputFrame.fftPoints;
    

    % Background Noise threshold estimation.

        %TPSW METHOD
            parametersTPSW = {};
            parametersTPSW.lengthSW = 51;
            parametersTPSW.gapSizeSW = 8;
            parametersTPSW.rejectionFactor = 4;
            parametersTPSW.deltaTPSW = 0.005; % THIS MUST BE IN dB.

            spectrumThresholdTPSW = PeakThreshold_TPSW(inputFrame,parametersTPSW,0);

        %SSE METHOD

            numberCoeffsSSE = 101;

            spectrumThreshold = PeakThreshold_SSE(inputFrame,numberCoeffsSSE,DEBUG);
        
    
    % Peak Detection

        %Finding all local maxima
            initialPeaks = NaN(1,totalFreqBins); 

            for freqCounter = 2:(totalFreqBins-1)
                if (powerSpectrumDB (freqCounter) > powerSpectrumDB (freqCounter-1) && powerSpectrumDB (freqCounter) > powerSpectrumDB (freqCounter+1))
                    
                initialPeaks(freqCounter) = powerSpectrumDB(freqCounter);
                end
            end

        %Eliminating everything below the threshold
            detectedPeaks = NaN(1,totalFreqBins);
            detectedPeakPositions = find(isfinite(initialPeaks));

            for freqCounter = detectedPeakPositions
                if initialPeaks(freqCounter) > spectrumThreshold(freqCounter)
                    detectedPeaks(freqCounter) = initialPeaks(freqCounter);
                end
            end

    % Enhancing the peaks and eliminating NaN elements to output arrays of size equal to the number of detected peaks.

        ENHANCEMENT = 1; 

        if ENHANCEMENT == 1
            detectedPeakPositions = find(isfinite(detectedPeaks));
            detectedPeakMatrix = PeakEnhancement(powerSpectrumDB,detectedPeaks,detectedPeakPositions,samplingRate,fftPoints);
            detectedFinalPeaks = detectedPeakMatrix(1,:);
            detectedPeakFrequencies = detectedPeakMatrix(2,:);
        else
            detectedPeakPositions = find(isfinite(detectedPeaks));
            detectedFinalPeaks = detectedPeaks(detectedPeakPositions);
            detectedPeakFrequencies = transpose(freqComponents(detectedPeakPositions));
        end


    if DEBUG == 1  

        figure
        hold on;
        plot(freqComponents,powerSpectrumDB,'G');
        plot(freqComponents,spectrumThreshold, 'R');
        plot(freqComponents,spectrumThresholdTPSW,'B');

        for freqCounter = detectedPeakPositions
            plot(freqComponents(freqCounter),detectedPeaks(freqCounter),'r*');
        end

        if ENHANCEMENT == 1

            fprintf('\nNumber of not-enhanced peaks: %i\n',length(detectedPeakPositions));
            fprintf('Number of enhanced peaks: %i power and %i frequency\n\n',length(detectedFinalPeaks),length(detectedPeakFrequencies));

            for freqCounter = 1:length(detectedPeakPositions)
                plot(detectedPeakFrequencies(freqCounter),detectedFinalPeaks(freqCounter),'k*');
            end

        end

        X = sprintf('Peak Detection of frame %i of %i',currentFrame,totalFrames);
            title(X);
            xlabel('Frequency (kHz)');
            ylabel('Power (dB)');
            legend ('Power Spectrum','SSE Threshold','TPSW Threshold','Detected Peaks(SSE)');

            hold off;
    
    end


end