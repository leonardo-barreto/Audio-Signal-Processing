function newTrack = createNewTrack(peakPower,peakFrequency,currentFrame)
    newTrack = {};
    newTrack.powerEvolution = peakPower; % Contains the power values of the track through its existence.
    newTrack.currentPower = peakPower;
    newTrack.frequencyEvolution = peakFrequency; % Contains the frequency values of the track through its existence.
    newTrack.startFrame = currentFrame; % Starting frame of the track (0 if non-existent)
    newTrack.finalFrame = 0; % Ending frame of the track (0 if non-existent or active)
    newTrack.length = 1; % length of the track (0 if non-existent)
    newTrack.hysteresis = 0; % Hysteresis counter
    newTrack.status = 'active'; % Track Status: 0 = active, 1 = inactive or 2 = asleep.
end 