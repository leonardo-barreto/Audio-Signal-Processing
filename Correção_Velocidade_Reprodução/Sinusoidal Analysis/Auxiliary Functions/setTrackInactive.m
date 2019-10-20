function updatedTrack = setTrackInactive(existingTrack,currentFrame,lastFrame)

    if (currentFrame == 1)
        error('You cannot update a track in the first frame, only create!' );
    elseif (existingTrack.status ~= 'asleep' && (currentFrame ~= lastFrame))
        error('Hmm... trying to deactivate a track that is not asleep before the final frame. Something is wrong.');
    end

    %Deactivacting a track that reached hysteresis limit or that didn't find a match in the last frame.

    existingTrack.hysteresis = NaN;
    existingTrack.status = 'inactive';

    updatedTrack = existingTrack;

end 