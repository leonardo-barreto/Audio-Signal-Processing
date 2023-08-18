function [existingTrack,deactivatedFlag] = treatUnmatchedTrack(existingTrack,maxHysteresis,currentFrame,lastFrame);

    if currentFrame == 1
        error('You cannot update a track in the first frame, only create!' );
    end

    deactivatedFlag = 0;

    switch existingTrack.status
        case 1 %Track is active and must go to sleep 
            existingTrack = setTrackAsleep(existingTrack,currentFrame);
        case 2 %Track is asleep
            if (existingTrack.hysteresis >= maxHysteresis) %Track is over hysteresis limit and must be deactivated.
                existingTrack = setTrackInactive(existingTrack,currentFrame,lastFrame);
                deactivatedFlag = 1;
            else %Track is asleep and can continue to be so.
                existingTrack = setTrackAsleep(existingTrack,currentFrame);
            end
        case 0 %Track is inactive
            error('Hmm... trying to treat a track that is inactive. Something is wrong.');
    end

end 