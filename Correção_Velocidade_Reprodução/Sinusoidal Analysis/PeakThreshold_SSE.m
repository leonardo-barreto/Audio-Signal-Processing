function spectrumFinalThreshold = PeakThreshold_SSE(inputFrame,numberCoeffsSSE,DEBUG);

    % This function takes a signal frame's power spectrum and executes a Spectrum Stochastic Estimating
    % process aiming at constructing a threshold for proeminent peak detection.
    %

    %Gathering frame data
    currentFrame = inputFrame.currentFrame;
    freqComponents = inputFrame.freqComponents;
    totalFreqBins = inputFrame.totalFreqBins;
    powerSpectrumDB = inputFrame.powerSpectrumDB;


    %This stage extends the spectrum by 20% of its size in order to avoid border effects during filtering.
    mirrorLength = floor(totalFreqBins/5);
    startMirror = flipud(powerSpectrumDB(1:mirrorLength));
    endMirror = flipud(powerSpectrumDB((totalFreqBins-mirrorLength)+1:totalFreqBins));
    powerSpectrumDB_extended = [startMirror;powerSpectrumDB;endMirror];
    extendedLength = length(powerSpectrumDB_extended); % length of extended spectrum.

    %Actual filtering
    spectrumS_filtered = movmean(powerSpectrumDB_extended,3);
    
    spectrumR_inverted = 1./spectrumS_filtered;

    spectrumR_inverted_filtered = movmean(spectrumR_inverted,numberCoeffsSSE);

    spectrumFinalThreshold = 1./spectrumR_inverted_filtered;

    %Discarding the mirrored edges and returning back to original size
    spectrumFinalThreshold = spectrumFinalThreshold((length(startMirror)+1):extendedLength-length(endMirror));

end