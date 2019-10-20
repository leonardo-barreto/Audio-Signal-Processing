function newTrack = setNewTrack(peakParameters,currentFrame)
    newTrack = {};
    newTrack.powerEvolution = peakParameters(1); % Contains the power values of the track through its existence.
    newTrack.currentPower = peakParameters(1);
    newTrack.frequencyEvolution = peakParameters(2); % Contains the frequency values of the track through its existence.
    newTrack.startFrame = currentFrame; % Starting frame of the track (0 if non-existent)
    newTrack.finalFrame = 0; % Ending frame of the track (0 if non-existent or active)
    newTrack.length = 1; % length of the track (0 if non-existent)
    newTrack.hysteresis = 0; % Hysteresis counter
    newTrack.status = 'active'; % Track Status: 0 = active, 1 = inactive or 2 = asleep.
end 