function detectedFinalPeaks = DetectSpectralPeaks(inputFrameSpectrum,totalFreqBins,samplingRate,windowSize,DEBUG)

    % This function aims to detect spectral peaks in a signal's frame, given its spectrum.
    %
    % Firstly, an initial all-peaks approach is taken.
    %
    % Secondly, a variable-length TPSW filtering builds a variable threshold.
    %
    % Finally, the frequency values are enhanced by use of the DFT1 method.
    %


    % Starts with an all-zero. This way, peaks have their correct bin number in the end.
    initialPeaks = zeros(1,1:totalFreqBins); 

    for freqCounter 2:(totalFreqBins-1)
        if inputFrameSpectrum (freqCounter) > inputFrameSpectrum (freqCounter-1) && 
           inputFrameSpectrum (freqCounter) > inputFrameSpectrum (freqCounter+1)
            
           initialPeaks(freqCounter) = inputFrameSpectrum(freqCounter);

        end
    end

    if DEBUG == 1
        %Plots all peaks.
        figure
        for freqCounter 1:totalFreqBins
            if (initialPeaks(freqCounter))
                stem(freqCounter,initialPeaks(freqCounter));
                hold on;
            else
                plot(freqCounter,inputFrameSpectrum(freqCounter));
                hold on;
            end
        end
        title('Initially Detected Peaks');
        xlabel('Frequency');
        ylabel('Power');
        hold off;
    end

    % TPSW starting point.
    



    % Frequency enhancement using the DFT1
    enhancedPeaks_DFT1 = EnhanceFrequency_DFT1(enhancedPeaks_TPSW,samplingRate);

    detectedFinalPeaks = enhancedPeaks_DFT1;
    
end