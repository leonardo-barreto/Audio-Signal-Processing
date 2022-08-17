function spectrumFinalThreshold = PeakThreshold_TPSW(inputFrame,parametersTPSW,DEBUG);

    % This function takes a signal frame's power spectrum and executes a Two-Pass Split Window
    % filtering process aiming to extract the most proeminent peaks of the spectrum, 
    % ignoring spurious or noisy components.

    %Gathering frame data
    currentFrame = inputFrame.currentFrame;
    freqComponents = inputFrame.freqComponents;
    totalFreqBins = inputFrame.totalFreqBins;
    powerSpectrum = inputFrame.powerSpectrum;

    %Gathering TPSW data
    lengthSW = parametersTPSW.lengthSW;
    gapSizeSW = parametersTPSW.gapSizeSW;
    rejectionFactor = parametersTPSW.rejectionFactor;
    deltaTPSW = parametersTPSW.deltaTPSW;


    %Building the Two-Pass Split Window
    if(gapSizeSW >= lengthSW)
        error('Invalid TPSW window. Length must be bigger than middle gap size.')
    end

    totalLengthSW = 2*lengthSW-1; % total length of the window
    middlePoint = floor((2*lengthSW-1)/2) + 1; %position of the center zero sample.
    nonZeroSamples = 2*(lengthSW-gapSizeSW); % number of non-zero samples

    splitWindow(middlePoint) = 0;

    for counter = 1:(lengthSW-1)
        if counter < gapSizeSW
            splitWindow(middlePoint+counter)=0;
            splitWindow(middlePoint-counter)=0;
        else
            splitWindow(middlePoint+counter)=1/nonZeroSamples;
            splitWindow(middlePoint-counter)=1/nonZeroSamples;
        end
    end


    %TPSW Filtering

    %This stage extends the spectrum by 20% of its size in order to avoid border effects during filtering.
    mirrorLength = floor(totalFreqBins/5);
    startMirror = flipud(powerSpectrum(1:mirrorLength));
    endMirror = flipud(powerSpectrum((totalFreqBins-mirrorLength)+1:totalFreqBins));
    powerSpectrum_extended = [startMirror;powerSpectrum;endMirror];

    %Actual filtering
    spectrumTPSW = conv(powerSpectrum_extended,splitWindow,'same'); %TPSW Filtering
    extendedLength = length(spectrumTPSW);


    %Substitution criterion
    for counter = 1:extendedLength
        if powerSpectrum_extended(counter) > rejectionFactor*spectrumTPSW(counter)
            spectrumSubstituted(counter) = spectrumTPSW(counter);
        else
            spectrumSubstituted(counter) = powerSpectrum_extended(counter);
        end
    end


    %Moving average filter
    spectrumFiltered = movmean(spectrumSubstituted,nonZeroSamples) + deltaTPSW;

    %Discarding the mirrored edges and returning back to original size (first two are for debug purposes)
    spectrumTPSW = spectrumTPSW((length(startMirror)+1):extendedLength-length(endMirror));
    spectrumSubstituted = spectrumSubstituted((length(startMirror)+1):extendedLength-length(endMirror));
    spectrumFinalThreshold = spectrumFiltered((length(startMirror)+1):extendedLength-length(endMirror));

    if DEBUG == 1

        figure;
        n = -lengthSW+1:1:lengthSW-1;
        stem(n,splitWindow,'filled');
        X = sprintf('Total length of split window: %i (N = %i and M = %i)',totalLengthSW,lengthSW,gapSizeSW);
        title(X);
    end   
    
end