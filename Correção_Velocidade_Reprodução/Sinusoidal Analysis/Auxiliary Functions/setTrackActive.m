function updatedTrack = setTrackActive(existingTrack,peakParameters,currentFrame)

    if (currentFrame == 1)
        error('You cannot update a track in the first frame, only create!' );
    end

    if (existingTrack.status == 2) %Bringing back a track that was asleep.

        hysteresisCount = existingTrack.hysteresis;
        
        existingTrack.powerEvolution(end+1:end+hysteresisCount) = existingTrack.powerEvolution(end);
        existingTrack.frequencyEvolution(end+1:end+hysteresisCount) = existingTrack.frequencyEvolution(end);
        existingTrack.hysteresis = 0;
        existingTrack.powerEvolution(end+1) = peakParameters(1);
        existingTrack.currentPower = peakParameters(1);
        existingTrack.frequencyEvolution(end+1) = peakParameters(2);
        existingTrack.currentFrequency = peakParameters(2);
        existingTrack.finalFrame = currentFrame;
        existingTrack.length = existingTrack.length + hysteresisCount + 1;
        existingTrack.status = 1;

    else %Updating a track that found a match

        if existingTrack.status ~= 1
            error('Hmm... trying to put a track that is not active or asleep into active. Something is wrong.');
        end

        existingTrack.powerEvolution(end+1) = peakParameters(1);
        existingTrack.currentPower = peakParameters(1);
        existingTrack.frequencyEvolution(end+1) = peakParameters(2);
        existingTrack.currentFrequency = peakParameters(2);
        existingTrack.finalFrame = currentFrame;
        existingTrack.length = existingTrack.length + 1;

    end

    updatedTrack = existingTrack;

end 