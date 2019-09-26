function y = TPSW_Filtering(inputFrameSpectrum,lengthSW,gapSizeSW,rejectionFactor);

    %Building the TPSW window

    if(gapSizeSW >= lengthSW)
        error('Invalid TPSW window. Length must be bigger than middle gap size.')
    end

    totalLength = lengthSW+ gapSizeSW+1; %overall length of TPSW
    middlePoint = floor(totalLength/2) + 1; %index of the zero in the middle.

    splitWindow(middlePoint) = 0;

    for counter = 1:(lengthSW-1)
        if counter < gapSizeSW
            splitWindow(middlePoint+counter)=0;
            splitWindow(middlePoint-counter)=0;
        else
            splitWindow(middlePoint+counter)=1;
            splitWindow(middlePoint-counter)=1;
        end
    end

    totalFreqBins = length(inputFrameSpectrum);

    TPSW_FilteredSequence = zeros(1,totalFreqBins);
    for counter = 1:totalFreqBins-(totalLength)
        sample = 0;
        for internal = 1:totalLength
            sample = sample + splitWindow(internal)*inputFrameSpectrum(counter+(internal-1));
        end
        TPSW_FilteredSequence(counter) = sample/totalFreqBins;
    end

    y = TPSW_FilteredSequence;
        
        



    
    
end