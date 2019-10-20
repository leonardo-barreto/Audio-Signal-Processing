function enhancedPeakMatrix = PeakEnhancement(powerSpectrumDB,peakPower,peakPositions,samplingRate,fftPoints)

    enhancedPeakMatrix = [];


    %Parabolic interpolation
        for peakIndex = 1:length(peakPositions);

            kp = peakPositions(peakIndex); %position of current peak (in original spectrum).

            A1 = powerSpectrumDB(kp-1); %Power of the sample to the left of the peak.
            A2 = peakPower(peakIndex); %Power of the current peak.
            A3 = powerSpectrumDB(kp+1); %Power of the sample to the right of the peak.

            d = (A1-A3)/2*(A1-2*A2+A3); %factor of adjustment to the peak position.

            newPeakFrequency = ((kp + d)*samplingRate)/fftPoints;
            newPeakPower = A2 - ((d*(A1-A3))/4);

            enhancedPeakMatrix(:,peakIndex) = [newPeakPower;newPeakFrequency];

        end

    %DFT1
        %enhancedFrequencies = (samplingRate/pi)*sin(pi*enhancedPeakMatrix(2,:)./samplingRate);
        %enhancedPeakMatrix(2,:) = enhancedFrequencies;


   

end