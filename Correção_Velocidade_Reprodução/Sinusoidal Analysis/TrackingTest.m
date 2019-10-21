function currentTracks = TrackingTest(signalFrame,currentTracks);


% ---------------------------------------- Sinusoidal Tracking -----------------------------------------------

%fprintf('\nSinusoidal tracking starting...\n');

% Array that holds all of the signal's tracks.
% This array starts with an empty track, and is extended as the signal demands new tracks, maintaining
% the information of tracks that ended. Later, another part of the algorithm will gather the starting
% and ending frames of each track and organize them.
if signalFrame.currentFrame == 1
    currentTracks = setNewTrack();
end

% Extending the signal frame built earlier to include some tracking parameters.
%signalFrame.currentFrame = 1; % Reset to signal start

%for frameCounter = 1:totalFrames

    %signalFrame.currentFrame = frameCounter;
    %signalFrame.peakMatrix = peakMatrix{frameCounter};

    %fprintf('\nI am at frame: %i and I have % tracks.\n',signalFrame.currentFrame,length(currentTracks));

    currentTracks = PartialTracking(signalFrame,currentTracks,1);

%end

