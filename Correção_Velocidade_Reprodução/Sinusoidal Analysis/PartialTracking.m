function currentTracks = PartialTracking(inputFrame,currentTracks,DEBUG)

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Defining parameters -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        maxTracksPerFrame = 100;

        freqTolerance = power(2,1/24); %quarter-tone (in Hz)
        maxHysteresis = 3; % in frames.
        minLength = 10; % in frames
        maxTrackFrequency = 20000; %in Hz
        minTrackPower = -50;%in dB

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Initial Processing -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|    

        %POWER in ROW 1, FREQUENCY in ROW 2
            powerRow = 1;
            freqRow = 2;

        % Gathering frame data
            currentFrame = inputFrame.currentFrame;
            peakMatrix = inputFrame.peakMatrix;
            lastFrame = inputFrame.totalFrames;


    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Gathering all frame peak information -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|


        %Discarding peaks above the max allowed frequency
            peakMatrix = sortcolumns(peakMatrix,freqRow,'ascend');

            allowedPeaks = peakMatrix(freqRow,:) < maxTrackFrequency;
            peakMatrix = peakMatrix(:,allowedPeaks);


        %Discarding peaks below minimum allowed power.
            peakMatrix = sortcolumns(peakMatrix,powerRow,'descend');
            allowedPeaks = peakMatrix(powerRow,:) > minTrackPower;
            peakMatrix = peakMatrix(:,allowedPeaks);

            if DEBUG == 1
                fprintf('%i th Frame information gathered.\n',currentFrame);
                fprintf('Peaks below max frequency and above minimum power: %i of %i\n',size(peakMatrix,2),size(inputFrame.peakMatrix,2));
                fprintf('Proceeding to update existing tracks...\n\n');
            end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Processing EXISTING Tracks -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        totalInactiveTracks = length(structArrayOperations(currentTracks,'status','==','inactive'));
        
        if (~isempty(currentTracks) && totalInactiveTracks ~= length(currentTracks)) % Will not attempt to process existing tracks if there aren't any.


            % -------------------------------- Matching tracks to peaks ------------------------------------
        
                activeIndexes = structArrayOperations(currentTracks,'status','==','active');  % Gathering active indexes
                asleepIndexes = structArrayOperations(currentTracks,'status','==','asleep'); %Gathering asleep indexes
                existingIndexes = [activeIndexes,asleepIndexes];  % Indexes of tracks to be continued: Active tracks are priority.

                if DEBUG == 1
                    fprintf('Existing (active or asleep) tracks: %i of %i.\n',length(existingIndexes),length(currentTracks));
                    fprintf('Proceeding to match tracks to peaks...\n\n');
                end

                %Actual matching
                unmatchedIndexes = []; %This will keep indexes of the tracks that didn't find a match.
                existingTrackFrequencies = getArrayFields(currentTracks,'currentFrequency',existingIndexes);
                trackIndex = 1;

                while (size(peakMatrix,2) > 0 && trackIndex <= length(existingIndexes))

                    freqTolerance = existingTrackFrequencies(trackIndex)*(freqTolerance-1);
                    %lowerBound = existingTrackFrequencies(trackIndex)/freqTolerance;

                    %biggerIndexes = find (peakMatrix(freqRow,:) >= lowerBound);
                    %smallerIndexes = find (peakMatrix(freqRow,:) <= upperBound);

                    %matchIndexes = intersect(biggerIndexes,smallerIndexes);

                    freqOffsets = abs(existingTrackFrequencies(trackIndex)-peakMatrix(freqRow,:));
                    [minOffset,matchIndex]= min(freqOffsets);

                    if (minOffset < freqTolerance)

                        if DEBUG == 1
                            fprintf('Track Match! This is track %i.\n',existingIndexes(trackIndex));
                        end
                        currentTracks(existingIndexes(trackIndex)) = setTrackActive(currentTracks(existingIndexes(trackIndex)),peakMatrix(:,matchIndex),currentFrame);
                        peakMatrix(:,matchIndex) = [];

                    else %no peaks are close enough.

                        if DEBUG == 1
                            fprintf('Track NOT Matched. This is track %i.\n',existingIndexes(trackIndex));
                        end

                        unmatchedIndexes(end+1) = existingIndexes(trackIndex);

                    end
                    trackIndex = trackIndex + 1;

                end

                if DEBUG == 1
                    fprintf('I have looked for tracks %i times.\n',trackIndex-1);
                end

                if DEBUG == 1
                    fprintf('Matching done.\n');
                    fprintf('Unmatched tracks: %i of %i previously existing.\n',length(unmatchedIndexes),length(existingIndexes));
                    fprintf('Proceeding to update unmatched tracks'' statuses...\n\n');
                end

            % ------------------------------------ Updating status of all unmatched tracks ---------------------------------------

                inactives = []; % This will index into existingIndexes later to remove inactive tracks.
                
                for trackIndex = 1:length(unmatchedIndexes)
                    if (currentTracks(unmatchedIndexes(trackIndex)).status == 1) %Track is active and must go to sleep. 
                        currentTracks(unmatchedIndexes(trackIndex)) = setTrackAsleep(currentTracks(unmatchedIndexes(trackIndex)),currentFrame);
                        if DEBUG == 2
                            fprintf('\nTrack %i is going to sleep.\n',unmatchedIndexes(trackIndex));
                        end
                    else

                        if (currentTracks(unmatchedIndexes(trackIndex)).status ~= 2)
                            error('Hmmm... unmatched track that is inactive? This should not happen.');
                        end

                        if (currentTracks(unmatchedIndexes(trackIndex)).hysteresis >= maxHysteresis) %Track is over hysteresis limit and must be deactivated.
                            currentTracks(unmatchedIndexes(trackIndex)) = setTrackInactive(currentTracks(unmatchedIndexes(trackIndex)),currentFrame,lastFrame);
                            if DEBUG == 2
                                fprintf('\nTrack %i is being deactivated.\n',unmatchedIndexes(trackIndex));
                            end
                            inactives(end+1) = trackIndex;
                        else %Track is asleep and can continue to be so.
                            currentTracks(unmatchedIndexes(trackIndex)) = setTrackAsleep(currentTracks(unmatchedIndexes(trackIndex)),currentFrame);
                            if DEBUG == 2
                                fprintf('\nTrack %i kept asleep.\n',unmatchedIndexes(trackIndex));
                            end
                        end
                    end 
                end

            % -------------------------------------- Updating track indexes ------------------------------------------

                existingIndexes(inactives) = [];
                activeIndexes = structArrayOperations(currentTracks,'status','==','active');

                if DEBUG == 1
                    asleepIndexes_debug = structArrayOperations(currentTracks,'status','==','asleep');
                    fprintf('Status updating done.\n');
                    fprintf('Deactivated tracks: %i\n',length(inactives));
                    fprintf('Asleep tracks: %i\n',length(asleepIndexes_debug));
                    fprintf('Active tracks: %i\n',length(activeIndexes));
                    fprintf('Existing tracks (asleep + active) after deactivations: %i\n',length(existingIndexes));
                    fprintf('Total tracks: %i\n',length(currentTracks));
                    if ((length(activeIndexes)+length(asleepIndexes_debug))~=length(existingIndexes))
                        error('Hmm... active + asleep is not equal to existing. Something is wrong and you dont know how to fix it...');
                    end
                    fprintf('Do all these make sense? asleep + active = existing is always necessary, so is existing + deactivated = total.\n')
                end
        
        else
            if DEBUG == 1
                fprintf('There are no existing tracks or all tracks are inactive. Therefore no tracks were updated.\n');
                fprintf('This is frame: %i.\n',currentFrame);
                %fprintf('If this is not frame 1, then either your tracks are coming up too short or there is something wrong.\n\n');
            end
        end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Allocating new tracks if allowed -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-||

        %Re-gathering total tracks
        activeIndexes = structArrayOperations(currentTracks,'status','==','active');  % Gathering active indexes
        totalActiveTracks = length(activeIndexes);
        totalTracks = length(currentTracks);

        if DEBUG == 1
            fprintf('This should allow for %i new tracks.\n',maxTracksPerFrame-totalActiveTracks);
            fprintf('There are %i allowed peaks.\n',size(peakMatrix,2));
        end

        peakIndex = 1;

        %New tracks should be allowed if: there is space for new tracks; max length of a new track is greater than the minimum allowed track length.
        while ((isempty(currentTracks)||(totalActiveTracks < maxTracksPerFrame))&&(lastFrame-currentFrame)>=(minLength-1))

            if (peakIndex > size(peakMatrix,2)) %ends the loop if there are no more peaks.
                break;
            end
            
            currentTracks(totalTracks + 1) = setNewTrack(peakMatrix(:,peakIndex),currentFrame);
            
            totalTracks = totalTracks + 1;
            totalActiveTracks = totalActiveTracks + 1;
            peakIndex = peakIndex + 1;

        end

        if DEBUG == 1
            fprintf('I have allocated %i tracks.\n',peakIndex-1);

    % ----------------------- Removing inactive tracks that ended too short -----------------------------

        inactiveIndexes = structArrayOperations(currentTracks,'status','==','inactive');
        if length(inactiveIndexes) > 0
            inactiveTrackLengths = getArrayFields(currentTracks,'length',inactiveIndexes);
            removeLogic = inactiveTrackLengths < minLength;

            deleteIndexes = inactiveIndexes(removeLogic);

            currentTracks(deleteIndexes) = [];

            if DEBUG == 1
                fprintf('Deleted %i of %i total inactive tracks.\n\n',length(deleteIndexes),length(inactiveIndexes));
            end
        else
            if DEBUG == 1
                fprintf('No tracks are inactive.\n\n');
            end
        end

        if currentFrame == lastFrame
            shortIndexes = structArrayOperations(currentTracks,'length','<',minLength);
            currentTracks(shortIndexes) = [];
        end

end

    

