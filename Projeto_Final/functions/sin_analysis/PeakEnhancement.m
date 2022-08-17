function enhancedPeakMatrix = PeakEnhancement(powerSpectrumDB,peakPower,peakPositions,samplingRate,fftPoints,DEBUG)

    enhancedPeakMatrix = [];


    %Parabolic interpolation
        if (~isempty(peakPositions))
            for peakIndex = 1:length(peakPositions);

                kp = peakPositions(peakIndex); %position of current peak (in original spectrum).

                A1 = powerSpectrumDB(kp-1); %Power of the sample to the left of the peak.
                A2 = powerSpectrumDB(kp); %Power of the current peak.
                A3 = powerSpectrumDB(kp+1); %Power of the sample to the right of the peak.

                d = (A1-A3)/2/(A1-2*A2+A3); %factor of adjustment to the peak position.

                newPeakFrequency = (((kp-1) + d)*samplingRate)/fftPoints;
                newPeakPower = A2 - ((d*(A1-A3))/4);

                enhancedPeakMatrix(:,peakIndex) = [newPeakPower;newPeakFrequency];
                %if DEBUG == 1
                %    fprintf('\n kp = %d. A1 = %d A2 = %d A3 = %d d = %d.\n',kp,A1,A2,A3,d);
                %    fprintf('newPower = %d. newFrequency = %d\n',newPeakPower,newPeakFrequency);
                %end

            end
        end

   

end