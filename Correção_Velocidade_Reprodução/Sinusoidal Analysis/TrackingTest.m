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
        signalFrame.peakMatrix = [detectedPeaksMatrix{frameCounter};detectedFrequenciesMatrix{frameCounter}];

        currentTracks = PartialTracking(signalFrame,currentTracks,1);
    end