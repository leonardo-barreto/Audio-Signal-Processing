function [frameArray] = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,spectrumSize,DEBUG)

    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - STFT
    %
    %   2nd - Frequency Enhancement
    %
    %   3rd - Peak Detection
    %
    %   4th - Construction of sinusoidal tracks
    %


    % -------------------------------------- STFT and spectrogram stage -------------------------------------------

        fftPoints = 2*spectrumSize; % Needs to account for double-sided spectrum.
        
        fprintf('\nSinusoidal Analysis started.\n Sampling Rate(Hz): %i\n Window: %s (size %i, overlap %i%%) \n Spectrum size: %i (%i FFT points)\n', samplingRate,windowType,windowSize,overlapPerc,spectrumSize,fftPoints);

        [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);
        
        powerMatrixDB = 10*log(powerMatrix); 

        totalFreqBins = length(freqComponents);
        totalFrames = length(timeInstants);

        % Building a signal frame as a peak detection entity
        signalFrame = {};
        signalFrame.totalFrames = totalFrames; %Total number of signal frames
        signalFrame.currentFrame = 1; %Current frame
        signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
        signalFrame.freqComponents = freqComponents; %frequency components vector
        signalFrame.timeInstants = timeInstants;
        signalFrame.samplingRate = samplingRate;
        signalFrame.fftPoints = fftPoints;

    % ---------------------------------------------- Peak Detection ------------------------------------------------

        fprintf('\nPeak Detection Starting...\n');

        detectedPeaksMatrix = {};
        detectedFrequenciesMatrix = {};

        if DEBUG == 1
            %Random frame chosen for DEBUG (temporary)
            DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);

            %DEBUG_FRAME = ra;

            for frameCounter = 1:totalFrames
                signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
                signalFrame.currentFrame = frameCounter;
                frameArray(frameCounter) = signalFrame;
                if DEBUG_FRAME == signalFrame.currentFrame
                    [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,1);
                else
                    frameArray(frameCounter) = signalFrame;
                    [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,0);
                end
            end

        else

            for frameCounter = 1:totalFrames
                signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
                signalFrame.currentFrame = frameCounter;
                frameArray(frameCounter) = signalFrame;
                [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,0);
            end

        end

        for index = 1:totalFrames
            frameArray(index).peakMatrix = [detectedPeaksMatrix{index};detectedFrequenciesMatrix{index}];
        end

    % ---------------------------------------- Sinusoidal Tracking -----------------------------------------------

        %fprintf('\nSinusoidal tracking starting...\n');

        % Extending the signal frame built earlier to include some tracking parameters.
            %signalFrame.currentFrame = 2; % Reset to signal start

        % Array that holds all of the signal's tracks.
        % This array starts with an empty track, and is extended as the signal demands new tracks, maintaining
        % the information of tracks that ended. Later, another part of the algorithm will gather the starting
        % and ending frames of each track and organize them.
            %if signalFrame.currentFrame == 1
                %currentTracks = setNewTrack();
            %end

        %Actual tracking
            %for frameCounter = 1:totalFrames
                %currentTracks = PartialTracking(frameArray(2),currentTracks,0);
            %end

end