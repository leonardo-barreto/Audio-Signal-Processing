function newTrack = setNewTrack(peakParameters,currentFrame)

    if nargin == 0 %Creates empty track
        newTrack.powerEvolution = []; % Contains the power values of the track through its existence.
        newTrack.currentPower = [];
        newTrack.frequencyEvolution = []; % Contains the frequency values of the track through its existence.
        newTrack.currentFrequency = [];
        newTrack.startFrame = 0; % Starting frame of the track (0 if non-existent)
        newTrack.finalFrame = 0; % Ending frame of the track (0 if non-existent or active)
        newTrack.length = 0; % length of the track (0 if non-existent)
        newTrack.hysteresis = 0; % Hysteresis counter
        newTrack.status = 0;
    else
        newTrack.powerEvolution = peakParameters(1); % Contains the power values of the track through its existence.
        newTrack.currentPower = peakParameters(1);
        newTrack.frequencyEvolution = peakParameters(2); % Contains the frequency values of the track through its existence.
        newTrack.currentFrequency = peakParameters(2);
        newTrack.startFrame = currentFrame; % Starting frame of the track (0 if non-existent)
        newTrack.finalFrame = currentFrame; % Ending frame of the track (0 if non-existent or active)
        newTrack.length = 1; % length of the track (0 if non-existent)
        newTrack.hysteresis = 0; % Hysteresis counter
        newTrack.status = 2; % Track Status: 0 = inactive, 1 = active, 2 = asleep.
    end

end 