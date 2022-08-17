function updatedTrack = setTrackAsleep(existingTrack,currentFrame)

    if (currentFrame == 1)
        error('You cannot update a track in the first frame, only create!' );
    end

    existingTrack.hysteresis = existingTrack.hysteresis + 1; % Incrementing hysteresis (counts for both active and asleep tracks).

    if (existingTrack.status ~= 2) % Putting an active track to sleep.
        
        if existingTrack.status ~= 1
            error('Hmm... trying to put a track that is not active or asleep to sleep. Something is wrong.');
        end

        existingTrack.status = 2;
        
    end

    updatedTrack = existingTrack;       

end 