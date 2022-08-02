function [detectedPeakMatrix,spectrumThreshold] = PeakDetection(inputFrame,samplingRate,fftPoints,DEBUG)

    % This function aims to detect spectral peaks in a signal's frame, given its spectrum.
    %
    % Firstly, a threshold estimation method is applied.
    %
    % Secondly, the peaks are detected based on the threshold.
    %
    % Finally, the peaks are submitted to a parabolic interpolation aiming to correct its frequency and power.
    %

    %Parameters
    peakProminence = 20; %in dB. This is for the findpeaks function. Controls how proeminent detected peaks must be.

    numberCoeffsSSE = 30;
    thresholdOffsetSSE = 0 ; %THIS MUST BE IN dB.

    hardThreshold = 80; %in DB
    freqThreshold = 5000; %in Hz

    % Gathering frame data
    
        powerSpectrum = inputFrame.powerSpectrum;
        freqComponents = inputFrame.freqComponents;
        totalFreqBins = inputFrame.totalFreqBins;
        totalFrames = inputFrame.totalFrames;
        currentFrame = inputFrame.currentFrame;

    

    % Background Noise threshold estimation.

        %TPSW METHOD
            parametersTPSW = {};
            parametersTPSW.lengthSW = 51;
            parametersTPSW.gapSizeSW = 8;
            parametersTPSW.rejectionFactor = 4;
            parametersTPSW.deltaTPSW = 5; % THIS MUST BE IN dB.

            spectrumThreshold = PeakThreshold_TPSW(inputFrame,parametersTPSW,0);

        %SSE METHOD
            spectrumThreshold_SSE = PeakThreshold_SSE(inputFrame,numberCoeffsSSE,DEBUG) + thresholdOffsetSSE;
        
    
    % Peak Detection

        %Finding all local maxima
            initialPeaks = NaN(1,totalFreqBins); 

            %for freqCounter = 2:(totalFreqBins-1)
            %    if (powerSpectrum (freqCounter) > powerSpectrum (freqCounter-1) && powerSpectrum (freqCounter) > powerSpectrum (freqCounter+1))
                    
            %    initialPeaks(freqCounter) = powerSpectrum(freqCounter);
            %    end
            %end

            [pks,locs] = findpeaks(powerSpectrum,'MinPeakProminence',peakProminence);

            initialPeaks(locs) = pks;


        %Eliminating everything below the threshold
            detectedPeaks = NaN(1,totalFreqBins);
            detectedPeakPositions = find(isfinite(initialPeaks));

            maxPeakPower = max(initialPeaks);

            for freqCounter = detectedPeakPositions
                if (initialPeaks(freqCounter) > spectrumThreshold(freqCounter) && (maxPeakPower - initialPeaks(freqCounter)) < hardThreshold && freqComponents(freqCounter) < freqThreshold)
                    detectedPeaks(freqCounter) = initialPeaks(freqCounter);
                end
            end

        %Eliminating NaN elements
            detectedPeakPositions = find(isfinite(detectedPeaks));
            detectedFinalPeaks = detectedPeaks(detectedPeakPositions);
            detectedPeakFrequencies = freqComponents(detectedPeakPositions);

    % Enhancing the peaks and eliminating NaN elements to output arrays of size equal to the number of detected peaks.

        ENHANCEMENT = 1; 

        if(~isempty(detectedPeakPositions))

            if ENHANCEMENT == 1
                detectedPeakMatrix = PeakEnhancement(powerSpectrum,detectedPeaks,detectedPeakPositions,samplingRate,fftPoints,DEBUG);
                detectedFinalPeaks = detectedPeakMatrix(1,:);
                detectedPeakFrequencies = detectedPeakMatrix(2,:);
            else
                detectedPeakFrequencies = freqComponents(detectedPeakPositions);
            end

        end

        detectedPeakMatrix = [detectedFinalPeaks;detectedPeakFrequencies];


end