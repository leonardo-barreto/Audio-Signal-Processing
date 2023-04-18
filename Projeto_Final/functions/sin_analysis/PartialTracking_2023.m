function currentTracks = PartialTracking_2023(inputFrame,totalFrames,currentTracks,DEBUG)

    DEBUG = 0;

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Defining parameters -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        maxTracksPerFrame = 100;

        freqTolerance = (power(2,1/24)-1); % about 3% (quarter-tone)
        maxHysteresis = 3;                 % in frames
        minTrackLength = 10;               % in frames
        maxTrackFrequency = 10000;         % in Hz
        minTrackPower = -60;               % in dB

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Initial Processing -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|    

        %POWER in ROW 1, FREQUENCY in ROW 2
            powerRow = 1;
            freqRow = 2;

        % Gathering frame data
            currentFrame = inputFrame.currentFrame;
            peakMatrix = inputFrame.peakMatrix;
            lastFrame = totalFrames;


    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Gathering all frame peak information -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        if (~isempty(peakMatrix))
            %Discarding peaks above the max allowed frequency
                peakMatrix = sortcolumns(peakMatrix,freqRow,'ascend');
                allowedPeaks = peakMatrix(freqRow,:) < maxTrackFrequency;
                peakMatrix = peakMatrix(:,allowedPeaks);


            %Discarding peaks below minimum allowed power.
                peakMatrix = sortcolumns(peakMatrix,powerRow,'descend');
                allowedPeaks = peakMatrix(powerRow,:) > minTrackPower;
                peakMatrix = peakMatrix(:,allowedPeaks);
        end

        totalPeaks = size(peakMatrix,2);

        if DEBUG == 1
            fprintf('Frame information gathered.\n');
            fprintf('Peaks below max frequency and above minimum power: %i of %i\n',size(peakMatrix,2),size(inputFrame.peakMatrix,2));
            fprintf('Proceeding to update existing tracks...\n\n');
        end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Processing EXISTING Tracks -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        totalInactiveTracks = length(structArrayOperations(currentTracks,'status','==','inactive'));
        
        if (~isempty(currentTracks) && totalInactiveTracks ~= length(currentTracks)) % Will not attempt to process existing tracks if there aren't any.


            % -------------------------------- Matching TRACKS to PEAKS ------------------------------------
        
            % Gathering track info
            matchedIndexes = [];
            unmatchedIndexes = [];
            deactivatedIndexes = []; % This will index into existingIndexes later to remove deactivated tracks.
            activeIndexes = structArrayOperations(currentTracks,'status','==','active');  % Gathering active indexes
            asleepIndexes = structArrayOperations(currentTracks,'status','==','asleep'); % Gathering asleep indexes
            availableIndexes = [asleepIndexes,activeIndexes];  % Indexes of tracks to be continued: ASLEEP tracks are priority.
            %existingIndexes = availableIndexes;

            if DEBUG == 1
                fprintf('Available (active or asleep) tracks: %i of %i.\n',length(availableIndexes),length(currentTracks));
                fprintf('Existing peaks: %i\n',totalPeaks);
                fprintf('Proceeding to match tracks to peaks...\n\n');
            end 

            for trackIndex = availableIndexes

                trackFrequency = currentTracks(trackIndex).frequencyEvolution(end);
                matchedPeaks = intersect(find(peakMatrix(freqRow,:) < trackFrequency*(1 + freqTolerance)),find(peakMatrix(freqRow,:) > trackFrequency*(1 - freqTolerance)));
                
                if (~isempty(matchedPeaks)) % Track matched to a peak
                    [~,matchIndex] = min(abs(peakMatrix(freqRow,matchedPeaks)-trackFrequency));
                    matchIndex = matchedPeaks(matchIndex);

                    if DEBUG == 1
                        fprintf('Peak %i matched to track %i.\n',matchIndex,trackIndex);
                    end

                    currentTracks(trackIndex) = setTrackActive(currentTracks(trackIndex),peakMatrix(:,matchIndex),currentFrame);
                    matchedIndexes(end+1) = trackIndex;
                    peakMatrix(:,matchIndex) = []; % Removes from selection a peak that was used

                else % No peaks are close enough
                    if (currentTracks(trackIndex).status == 0)
                        error('Hmmm... trying to process track that is inactive? This should not happen.');
                    end
                    if DEBUG == 1
                        fprintf('Track %i NOT matched to any peak.\n',trackIndex);
                    end
                    unmatchedIndexes(end+1) = trackIndex
                    switch currentTracks(trackIndex).status
                        case 1 %Track is active and must go to sleep 
                            currentTracks(trackIndex) = setTrackAsleep(currentTracks(trackIndex),currentFrame);
                        case 2 %Track is asleep
                            if (currentTracks(trackIndex).hysteresis >= maxHysteresis) %Track is over hysteresis limit and must be deactivated.
                                currentTracks(trackIndex) = setTrackInactive(currentTracks(trackIndex),currentFrame,lastFrame);
                                deactivatedIndexes(end+1) = trackIndex;
                            else %Track is asleep and can continue to be so.
                                currentTracks(trackIndex) = setTrackAsleep(currentTracks(trackIndex),currentFrame);
                            end
                    end
                end 
            end
            
            if DEBUG == 1
                fprintf('Matching done.\n');
                fprintf('Matched tracks: %i of %i previously existing.\n',length(matchedIndexes),length(availableIndexes));
                fprintf('Unmatched tracks: %i of %i previously existing.\n',length(unmatchedIndexes),length(availableIndexes));
                fprintf('%i peaks remaining of %i.\n\n',size(peakMatrix,2),totalPeaks);
            end

                
            % -------------------------------------- Updating track indexes ------------------------------------------

                availableIndexes(deactivatedIndexes == availableIndexes) = [];
                existingIndexes = availableIndexes;
                activeIndexes = structArrayOperations(currentTracks,'status','==','active');

                if DEBUG == 1
                    asleepIndexes_debug = structArrayOperations(currentTracks,'status','==','asleep');
                    fprintf('Status updating done.\n');
                    fprintf('Deactivated tracks: %i\n',length(deactivatedIndexes));
                    fprintf('Asleep tracks: %i\n',length(asleepIndexes_debug));
                    fprintf('Active tracks: %i\n',length(activeIndexes));
                    fprintf('Existing tracks (asleep + active): %i\n',length(existingIndexes));
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
    % PAREI AQUI
    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Allocating new tracks if allowed -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-||

        %Re-gathering total tracks
        activeIndexes = structArrayOperations(currentTracks,'status','==','active');  % Gathering active indexes
        totalActiveTracks = length(activeIndexes);
        totalTracks = length(currentTracks);

        if DEBUG == 1
            fprintf('\nThis should allow for %i new tracks.\n',maxTracksPerFrame-totalActiveTracks);
        end

        peakIndex = 1;

        %New tracks should be allowed if: 
        %   - there is space for new tracks; 
        %   - max length of a new track is greater than the minimum allowed track length.
        
        while (((isempty(currentTracks)||(totalActiveTracks < maxTracksPerFrame))))

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
        end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Removing inactive tracks that ended too short -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        if currentFrame == lastFrame
            for index = 1:length(currentTracks)
                currentTracks(index) = setTrackInactive(currentTracks(index),currentFrame,lastFrame);
            end
        end

        inactiveIndexes = structArrayOperations(currentTracks,'status','==','inactive');

        if (~isempty(inactiveIndexes))

            inactiveTrackLengths = getArrayFields(currentTracks,'length',inactiveIndexes);
            removeLogic = inactiveTrackLengths < minTrackLength;

            deleteIndexes = inactiveIndexes(removeLogic);

            currentTracks(deleteIndexes) = [];

            if DEBUG == 1
                fprintf('Deleted %i of %i total inactive tracks.\n\n',length(deleteIndexes),length(inactiveIndexes));
            end
        else
            if DEBUG == 1
                fprintf('No tracks are inactive in frame %i.\n',currentFrame);
            end
        end

end

    

