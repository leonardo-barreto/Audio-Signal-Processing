function [y,x] = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,DEBUG)

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

    fprintf('\nSinusoidal Analysis started.\n Sampling Rate(Hz): %i\n Window: %s (size %i, overlap %i%%) \n FFT Points: %i\n', samplingRate,windowType,windowSize,overlapPerc,fftPoints);

    [spectrgMatrix,freqComponents_normalized,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);
    
    powerMatrixDB = 10*log(powerMatrix);
    freqComponents = (samplingRate/2*pi).*freqComponents_normalized;

    totalFreqBins = length(freqComponents);
    totalFrames = length(timeInstants);

    % Building a signal frame as a peak detection entity
    signalFrame = {};
    signalFrame.totalFrames = totalFrames; %Total number of signal frames
    signalFrame.currentFrame = 1; %Current frame
    signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
    signalFrame.freqComponents = freqComponents; %frequency components vector

    % ---------------------------------------------- Peak Detection ------------------------------------------------

    fprintf('\nPeak Detection Started.\n');

    detectedPeaksMatrix = {};
    detectedFrequenciesMatrix = {};

    if DEBUG == 1
        %Random frame chosen for DEBUG (temporary)
        DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
            signalFrame.currentFrame = frameCounter;
            if DEBUG_FRAME == signalFrame.currentFrame
                [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,1);
            else
                [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,0);
            end
        end

    else

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
            signalFrame.currentFrame = frameCounter;
            [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,0);
        end

    end

    % ---------------------------------------- Sinusoidal Tracking -----------------------------------------------

    %magnitudeMatrixDB = 10*log(abs(spectrgMatrix));

    % Building what is a template sinusoidal track
    templateTrack = {};
    templateTrack.powerEvolution = []; % Contains the power values of the track through its existence.
    templateTrack.currentPower = [];
    templateTrack.frequencyEvolution = []; % Contains the frequency values of the track through its existence.
    templateTrack.startFrame = 0; % Starting frame of the track (0 if non-existent)
    templateTrack.finalFrame = 0; % Ending frame of the track (0 if non-existent or active)
    templateTrack.length = 0; % length of the track (0 if non-existent)
    templateTrack.hysteresis = 0; % Hysteresis counter
    templateTrack.status = 'inactive'; % Track Status: 0 = active, 1 = inactive or 2 = asleep.
    
    % Array that holds all of the signal's tracks.
    % This array starts off at one track, and is extended as the signal demands new tracks, maintaining
    % the information of tracks that ended. Later, another part of the algorithm will gather the starting
    % and ending frames of each track and organize them.
    currentTracks = {};

    % Extending the signal frame built earlier to include some tracking parameters.
    signalFrame.currentFrame = 1; % Reset to signal start
    signalFrame.totalTracks = 0; % Total number of tracks in the frame.
    
    for frameCounter = 1:totalFrames

        signalFrame.currentFrame = frameCounter;
        signalFrame.spectrumPeaks = detectedPeaksMatrix{frameCounter};
        signalFrame.peakFrequencies = detectedFrequenciesMatrix{frameCounter};

        currentTracks = PartialTracking(signalFrame,currentTracks,1);
    end

    
    %TEMPORARY
    y = detectedPeaksMatrix;
    x = detectedPositionMatrix;

end