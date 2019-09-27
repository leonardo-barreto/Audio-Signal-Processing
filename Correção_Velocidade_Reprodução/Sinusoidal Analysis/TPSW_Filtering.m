function filteredSpectrum_finalThreshold = TPSW_Filtering(inputFrame,lengthSW,gapSizeSW,rejectionFactor,deltaTPSW,DEBUG,DEBUG_FRAME);

    % This function takes a signal frame's power spectrum and executes a Two-Pass Split Window
    % filtering process aiming to extract the most proeminent peaks of the spectrum, 
    % ignoring spurious or noisy components.

    %Gathering frame data
    currentFrame = inputFrame.currentFrame;
    freqComponents = inputFrame.freqComponents;
    totalFreqBins = inputFrame.totalFreqBins;
    powerSpectrumDB = inputFrame.powerSpectrumDB;


    %Building the TPSW window
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
    filteredSpectrum_TPSW = conv(powerSpectrumDB,splitWindow,'same'); %TPSW Filtering

    %Substitution criterion
    for counter = 1:totalFreqBins
        if powerSpectrumDB(counter) > rejectionFactor*filteredSpectrum_TPSW(counter)
            filteredSpectrum_substituted(counter) = filteredSpectrum_TPSW(counter);
        else
            filteredSpectrum_substituted(counter) = powerSpectrumDB(counter);
        end
    end

    %Moving average filter
    filteredSpectrum_finalThreshold = movmean(filteredSpectrum_substituted,nonZeroSamples) + deltaTPSW;

    if (DEBUG == 1 && DEBUG_FRAME == currentFrame)

        figure;
        n = -lengthSW+1:1:lengthSW-1;
        stem(n,splitWindow,'filled');
        X = sprintf('Total length of split window: %i (N = %i and M = %i)',totalLengthSW,lengthSW,gapSizeSW);
        title(X);

        figure;
        hold on;
        plot(freqComponents,powerSpectrumDB,'B');
        plot(freqComponents,filteredSpectrum_TPSW,'R');
        plot(freqComponents,filteredSpectrum_substituted,'G');
        plot(freqComponents,filteredSpectrum_finalThreshold,'K');
        X = sprintf ('TPSW Process of frame %i of %i',currentFrame,inputFrame.totalFrames);
        title(X);
        legend('Power Spectrum(DB)','TPSW-Filtered Spectrum','Substitution Results','Final Threshold');
        hold off;

    end   
    
end