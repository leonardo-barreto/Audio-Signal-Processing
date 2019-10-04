function [detectedFinalPeaks,detectedPeakPositions] = DetectSpectralPeaks(inputFrame,DEBUG)

    % This function aims to detect spectral peaks in a signal's frame, given its spectrum.
    %
    % Firstly, an initial all-peaks approach is taken.
    %
    % Secondly, a variable-length TPSW filtering builds a variable threshold.
    %
    % Finally, the frequency values are enhanced by use of the DFT1 method.
    %

    powerSpectrumDB = inputFrame.powerSpectrumDB;
    freqComponents = inputFrame.freqComponents;
    totalFreqBins = inputFrame.totalFreqBins;

    % Background Noise threshold estimation.

    % TPSW METHOD
    %parametersTPSW = {};
    %parametersTPSW.lengthSW = 50;
    %parametersTPSW.gapSizeSW = 5;
    %parametersTPSW.rejectionFactor = ;
    %parametersTPSW.deltaTPSW = 10; % THIS MUST BE IN dB.

    %spectrumThreshold = PeakThreshold_TPSW(inputFrame,parametersTPSW,DEBUG);

    %SSE METHOD

    numberCoeffsSSE = 50;

    spectrumThreshold = PeakThreshold_SSE(inputFrame,numberCoeffsSSE,DEBUG);
    
    % All-peaks detection (all local maxima)
    initialPeaks = zeros(1,totalFreqBins); 

    for freqCounter = 2:(totalFreqBins-1)
        if (powerSpectrumDB (freqCounter) > powerSpectrumDB (freqCounter-1) && powerSpectrumDB (freqCounter) > powerSpectrumDB (freqCounter+1))
            
           initialPeaks(freqCounter) = powerSpectrumDB(freqCounter);
        end
    end

    % Final peak detection

    detectedFinalPeaks = zeros(1,totalFreqBins);
    detectedPeakPositions = find(initialPeaks);

    for freqCounter = detectedPeakPositions
        if initialPeaks(freqCounter) > spectrumThreshold(freqCounter)
            detectedFinalPeaks(freqCounter) = initialPeaks(freqCounter);
        end
    end

end