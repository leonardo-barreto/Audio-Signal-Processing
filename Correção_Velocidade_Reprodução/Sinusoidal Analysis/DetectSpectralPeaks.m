function [detectedFinalPeaks,detectedPeakFrequencies] = DetectSpectralPeaks(inputFrame,samplingRate,fftPoints,DEBUG)

    % This function aims to detect spectral peaks in a signal's frame, given its spectrum.
    %
    % Firstly, a threshold estimation method is applied.
    %
    % Secondly, the peaks are detected based on the threshold.
    %
    % Finally, the peaks are submitted to a parabolic interpolation aiming to correct its frequency and power.
    %

    %Parameters
    peakProminence = 40; %in dB.

    numberCoeffsSSE = 20;
    thresholdOffsetSSE = 0 ; %THIS MUST BE IN dB.

    hardThreshold = 80; %in DB

    % Gathering frame data
    
        powerSpectrumDB = inputFrame.powerSpectrumDB;
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
            parametersTPSW.deltaTPSW = 30; % THIS MUST BE IN dB.

            spectrumThreshold_TPSW = PeakThreshold_TPSW(inputFrame,parametersTPSW,0);

        %SSE METHOD
            spectrumThreshold = PeakThreshold_SSE(inputFrame,numberCoeffsSSE,DEBUG) + thresholdOffsetSSE;
        
    
    % Peak Detection

        %Finding all local maxima
            initialPeaks = NaN(1,totalFreqBins); 

            %for freqCounter = 2:(totalFreqBins-1)
            %    if (powerSpectrumDB (freqCounter) > powerSpectrumDB (freqCounter-1) && powerSpectrumDB (freqCounter) > powerSpectrumDB (freqCounter+1))
                    
            %    initialPeaks(freqCounter) = powerSpectrumDB(freqCounter);
            %    end
            %end

            [pks,locs] = findpeaks(powerSpectrumDB,'MinPeakProminence',peakProminence);

            initialPeaks(locs) = pks;


        %Eliminating everything below the threshold
            detectedPeaks = NaN(1,totalFreqBins);
            detectedPeakPositions = find(isfinite(initialPeaks));

            maxPeakPower = max(initialPeaks);

            for freqCounter = detectedPeakPositions
                if (initialPeaks(freqCounter) > spectrumThreshold(freqCounter) && (maxPeakPower - initialPeaks(freqCounter)) < hardThreshold)
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
                detectedPeakMatrix = PeakEnhancement(powerSpectrumDB,detectedPeaks,detectedPeakPositions,samplingRate,fftPoints,DEBUG);
                detectedFinalPeaks = detectedPeakMatrix(1,:);
                detectedPeakFrequencies = detectedPeakMatrix(2,:);
            else
                detectedPeakFrequencies = freqComponents(detectedPeakPositions);
            end

        end


    if DEBUG == 1  
        figure
        hold on;
        plot(freqComponents./1000,powerSpectrumDB,'G','LineWidth',2);
        %plot(freqComponents./1000,spectrumThresholdSSE, 'R');
        plot(freqComponents./1000,spectrumThreshold,'B','LineWidth',2);

        %for freqCounter = detectedPeakPositions
        %    plot(freqComponents(freqCounter)/1000,detectedPeaks(freqCounter),'r*');
        %end

        if ENHANCEMENT == 1

            for freqCounter = 1:length(detectedPeakPositions)
                plot(detectedPeakFrequencies(freqCounter)/1000,detectedFinalPeaks(freqCounter),'r.','markersize',30);
            end

        end

        X = sprintf('Deteccao de picos do quadro %i de %i',currentFrame,totalFrames);
            title(X,'FontSize', 30);
            xlabel('Frequencia (kHz)','FontSize', 30);
            ylabel('Potencia (dB)','FontSize', 30);
            legend ('Espectro','Limiar SSE','Picos detectados');
            set(gca,'FontSize', 30);
            xlim([0 4.4])
            hold off;
    
    end


end