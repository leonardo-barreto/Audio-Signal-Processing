function spectrumFinalThreshold = PeakThreshold_SSE(inputFrame,numberCoeffsSSE,DEBUG);

    % This function takes a signal frame's power spectrum and executes a Spectrum Stochastic Estimating
    % process aiming at constructing a threshold for proeminent peak detection.
    %

    %Gathering frame data
    currentFrame = inputFrame.currentFrame;
    freqComponents = inputFrame.freqComponents;
    totalFreqBins = inputFrame.totalFreqBins;
    powerSpectrumDB = inputFrame.powerSpectrum;

    powerSpectrum = powerSpectrumDB/10;
    powerSpectrum = power(10,powerSpectrum);


    %This stage extends the spectrum by 20% of its size in order to avoid border effects during filtering.
    mirrorLength = floor(totalFreqBins/5);
    startMirror = flipud(powerSpectrum(1:mirrorLength));
    endMirror = flipud(powerSpectrum((totalFreqBins-mirrorLength)+1:totalFreqBins));
    powerSpectrum_extended = [startMirror;powerSpectrum;endMirror];
    extendedLength = length(powerSpectrum_extended); % length of extended spectrum.

    %Actual filtering
    spectrumS_filtered = movmean(powerSpectrum_extended,3);
    
    spectrumR_inverted = 1./spectrumS_filtered;

    spectrumR_inverted_filtered = movmean(spectrumR_inverted,numberCoeffsSSE);

    spectrumFinalThreshold = 1./spectrumR_inverted_filtered;

    %Discarding the mirrored edges and returning back to original size
    spectrumFinalThreshold = spectrumFinalThreshold((length(startMirror)+1):extendedLength-length(endMirror));

    spectrumFinalThreshold = 10*log10(spectrumFinalThreshold);

    %maxPower = max(powerSpectrumDB);
    %powerThreshold = maxPower - 80;

    %for index = 1:length(spectrumFinalThreshold)
        %correctedThreshold(index) = max(powerThreshold,spectrumFinalThreshold(index));
    %end

    %spectrumFinalThreshold = correctedThreshold;

end