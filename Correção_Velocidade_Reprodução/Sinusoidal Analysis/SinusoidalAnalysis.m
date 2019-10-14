function y = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,DEBUG)

    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - STFT
    %
    %   2nd - Frequency Enhancement
    %
    %   3rd - Peak Detection
    %
    %   4th - Tracking of sinusoidal components
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

    % ------------------------------------- Peak Detection ------------------------------------

    fprintf('\nPeak Detection Started.\n');

    detectedPeaksMatrix = zeros(totalFreqBins,totalFrames);

    if DEBUG == 1
        %Random frame chosen for DEBUG (temporary)
        DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            if DEBUG_FRAME == signalFrame.currentFrame
                detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,1);
            else
                detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,0);
            end
        end

    else

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,0);
        end

    end

    %TEMPORARY
    y = detectedPeaksMatrix;

    % ---------------------------- Sinusoidal Tracking ------------------------------------

    MAXTRACKSPERFRAME = 100; % Maximum of tracks per frame.
    magnitudeMatrixDB = 10*log(abs(spectrgMatrix));

    % Building what is a template sinusoidal track
    templateTrack = {};
    templateTrack.amplitudeEvolution = zeros(totalFrames); % Contains the amplitude values of the track through its existence.
    templateTrack.frequencyEvolution = zeros(totalFrames); % Contains the frequency values of the track through its existence.
    templateTrack.startFrame = 0; % Starting frame of the track (0 if non-existent)
    templateTrack.finalFrame = 0; % Ending frame of the track (0 if non-existent or active)
    templateTrack.length = 0; % length of the track (0 if non-existent)
    templateTrack.hysteresis = 0; % Hysteresis counter
    templateTrack.status = 'inactive'; % Track Status: inactive, active or asleep.

    %Creating one track
    singleTrack = templateTrack;

    % Array that holds all of the signal's tracks
    currentTracks(1) = singleTrack;

    % Extending the signal frame built earlier to include some tracking parameters.
    signalFrame.currentFrame = 1; % Reset to signal start
    signalFrame.totalTracks = 0; % Total number of tracks in the frame.
    
    for frameCounter = 1:totalFrames
        signalFrame.magnitudeSpectrumDB = (magnitudeMatrixDB(:,frameCounter));
        signalFrame.currentFrame = frameCounter;
        signalFrame.totalTracks = 0;
        currentTracks = PartialTracking(signalFrame,currentTracks,MAXTRACKSPERFRAME,1);
    end

    
    

end