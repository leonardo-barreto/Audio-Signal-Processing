function updatedTrack = updateTrack(existingTrack,peakParameters,currentFrame,newStatus)

    if (currentFrame == 1)
        error('You cannot update a track in the first frame, only create.' );
    end

    switch nargin

        case 3 % Track status will NOT be changed (Therefore, only active or asleep tracks go through this stage)

            existingTrack.powerEvolution(end+1) = peakParameters(1);
            existingTrack.currentPower = peakParameters(1);
            existingTrack.frequencyEvolution(end+1) = peakParameters(2);
            existingTrack.length = existingTrack.length + 1;

            if (newStatus == 'asleep')

            if (existingTrack.status == 'asleep')
                existingTrack.hysteresis = existingTrack.hysteresis + 1;
                existingTrack.powerEvolution(end+1) = powerEvolution(end);
                existingTrack.frequencyEvolution(end+1) = frequencyEvolution(end);
                existingTrack.length = existingTrack.length + 1;
            end

        case 4 % Track status WILL be changed

            existingTrack.hysteresis = existingTrack.hysteresis + 1;

            if (newStatus == 'asleep')
                existingTrack.powerEvolution(end+1) = peakParameters(1);
                existingTrack.currentPower = peakParameters(1);
                existingTrack.frequencyEvolution(end+1) = peakParameters2);
                existingTrack.length = existingTrack.length + 1;
                existingTrack.status = 'asleep';
            elseif (newStatus == 'inactive')
                existingTrack.finalFrame = (currentFrame-1);
                existingTrack.status = 'inactive';
            end

        otherwise
            error('Invalid number of arguments (min 3 max 4)');
    end

    updatedTrack = existingTrack;
            
end 