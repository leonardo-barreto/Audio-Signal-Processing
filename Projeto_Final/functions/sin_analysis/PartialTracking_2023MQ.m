function currentTracks = PartialTracking_2023MQ(inputFrame,TFParams,currentTracks,backwardsFlag)

    DEBUG = 0;

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Defining parameters -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        maxTracksPerFrame = 100;

        freqTolerance = (power(2,1/24)-1); % about 3% (quarter-tone)
        powerTolerance = 3;                % in dB
        maxHysteresis = 10;                % in frames
        minTrackLength = 20;               % in frames
        maxTrackFrequency = 5000;          % in Hz
        minTrackPower = -60;               % in dB

        %POWER in ROW 1, FREQUENCY in ROW 2
            powerRow = 1;
            freqRow = 2;
        % Gathering frame data
            currentFrame = inputFrame.currentFrame;
            peakMatrix = inputFrame.peakMatrix;
            totalFrames = length(TFParams.timeInstants);

            if backwardsFlag == 1
                currentFrame = totalFrames - (currentFrame - 1);
            end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Peak pre-processing -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        if (~isempty(peakMatrix))
            %Discarding peaks above the max allowed frequency and below minimum allowed power
            allowedPeaks = peakMatrix(freqRow,:) < maxTrackFrequency;
            peakMatrix = peakMatrix(:,allowedPeaks);
            allowedPeaks = peakMatrix(powerRow,:) > minTrackPower;
            peakMatrix = peakMatrix(:,allowedPeaks); 
            %Ordering peaks in frequency
            peakMatrix = sortcolumns(peakMatrix,freqRow,'ascend');
        end
        totalPeaks = size(peakMatrix,2);
        if DEBUG == 1
            fprintf('Frame information gathered.\n');
            fprintf('No. of peaks below max frequency and above minimum power: %i of %i\n',size(peakMatrix,2),size(inputFrame.peakMatrix,2));
            fprintf('Proceeding to update existing tracks...\n\n');
        end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Track pre-processing  -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        % Gathering track info
        inactiveIdxs = structArrayOperations(currentTracks,'status','==','inactive');
        inactiveTracks = currentTracks(inactiveIdxs);
        availableTracks = currentTracks(~ismember(1:length(currentTracks),inactiveIdxs)); % Tracks available for the matching process
        if DEBUG >= 1
            prevAvailableTracks = length(availableTracks);
            fprintf('\nAvailable (active or asleep) tracks: %i of %i.\n',length(availableTracks),length(currentTracks));
            fprintf('Existing peaks: %i\n',totalPeaks);
            fprintf('Proceeding to match tracks to peaks...\n\n');
        end 
        
    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Processing existing tracks -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
        
        % Won't attempt to process existing tracks if there aren't any.
        if (~isempty(availableTracks))
            % Creating necessary arrays
            matchedFlag = false(length(availableTracks),1);
            unmatchedFlag = false(length(availableTracks),1);
            deactivatedFlag = false(length(availableTracks),1);
            deleteIdxs = [];
            % Ordering tracks in frequency
            availableTracks = sortStruct(availableTracks, 'currentFrequency');
            for trackIdx = 1:length(availableTracks)
                trackFrequency = availableTracks(trackIdx).frequencyEvolution(end);
                matchedPeakIdxs = intersect(find(peakMatrix(freqRow,:) < trackFrequency*(1 + freqTolerance)),find(peakMatrix(freqRow,:) > trackFrequency*(1 - freqTolerance)));
                %Power matching criterion (experimental)
                %matchedPeaksPower = intersect(find(peakMatrix(powerRow,:) < availableTracks(trackIdx).powerEvolution(end) + powerTolerance),find(peakMatrix(powerRow,:) > availableTracks(trackIdx).powerEvolution(end) - powerTolerance));
                %matchedPeakIdxs = intersect(matchedPeakIdxs,matchedPeaksPower);

                if (~isempty(matchedPeakIdxs)) % 1) Track matched to a peak
                    [freqGap,peakCandidateIdx] = min(abs(peakMatrix(freqRow,matchedPeakIdxs)-trackFrequency));
                    peakCandidateIdx = matchedPeakIdxs(peakCandidateIdx);
                    if DEBUG == 1
                        fprintf('Peak %i matched to track %i.\n',peakCandidateIdx,trackIdx);
                    end
                    % 2) Check if there is better match to the peak
                    betterTrackIdx = find(abs(peakMatrix(freqRow,peakCandidateIdx)-getArrayFields(availableTracks(trackIdx+1:end),'currentFrequency')) < freqGap);
                    if (isempty(betterTrackIdx)) % current track is best match
                        availableTracks(trackIdx) = setTrackActive(availableTracks(trackIdx),peakMatrix(:,peakCandidateIdx),currentFrame);
                        matchedFlag(trackIdx) = 1;
                        peakMatrix(:,peakCandidateIdx) = [];
                        if DEBUG == 1
                            fprintf('Peak %i definite match to track %i.\n',peakCandidateIdx,trackIdx);
                        end  
                    else % 2.1) Better match for peak found
                        if DEBUG == 1
                            fprintf('Better match for peak %i found.\n',peakCandidateIdx,betterTrackIdx);
                        end
                        if (numel(matchedPeakIdxs) > 1 && peakCandidateIdx > 1) % 2.1.a) Track can still be matched to adjacent peak
                            availableTracks(trackIdx) = setTrackActive(availableTracks(trackIdx),peakMatrix(:,peakCandidateIdx-1),currentFrame);
                            matchedFlag(trackIdx) = 1;
                            peakMatrix(:,peakCandidateIdx-1) = [];
                            if DEBUG == 1
                                fprintf('Track %i matched definitely to peak %i instead.\n',trackIdx,peakCandidateIdx-1);
                            end
                        else % 2.1.b) Track can't be matched
                            unmatchedFlag(trackIdx) = 1;
                            [availableTracks(trackIdx),deactivatedFlag(trackIdx)] = treatUnmatchedTrack(availableTracks(trackIdx),maxHysteresis,currentFrame,totalFrames);
                            if DEBUG == 1
                                fprintf('No definite match possible for track %i.\n',trackIdx);
                            end
                        end 
                    end 
                else % 1) No matches found
                    unmatchedFlag(trackIdx) = 1;
                    [availableTracks(trackIdx),deactivatedFlag(trackIdx)] = treatUnmatchedTrack(availableTracks(trackIdx),maxHysteresis,currentFrame,totalFrames);
                    if DEBUG == 1
                        fprintf('Track %i NOT matched to any peak.\n',trackIdx);
                    end
                end 
            end
            % Deleting deactivated tracks that ended too short
            deactivatedIdxs = find(deactivatedFlag);
            if (~isempty(deactivatedIdxs))
                deleteIdxs = deactivatedIdxs(getArrayFields(availableTracks,'length',deactivatedIdxs) < minTrackLength);
                % Deleting tracks
                availableTracks(deleteIdxs) = [];
                % Removing deleted tracks' positions from flag vectors
                matchedFlag(deleteIdxs) = [];
                unmatchedFlag(deleteIdxs) = [];
                deactivatedFlag(deleteIdxs) = [];
                if DEBUG == 1
                    fprintf('Deleted %i of %i total inactive tracks.\n\n',length(deleteIdxs),length(inactiveIdxs));
                end
            else
                if DEBUG == 1
                    fprintf('No tracks are inactive in frame %i.\n',currentFrame);
                end
            end
            % Reorganizing available tracks only with active/asleep and separating inactives
            deactivatedIdxs = find(deactivatedFlag);
            inactiveTracks = [inactiveTracks availableTracks(deactivatedIdxs)];
            availableTracks = availableTracks(~ismember(1:length(availableTracks),deactivatedIdxs));
            % Status print (DEBUG)
            if DEBUG >= 1
                matchedIdxs = find(matchedFlag);
                unmatchedIdxs = find(unmatchedFlag);
                
                fprintf('\nMatching done.\n');
                fprintf('Matched tracks: %i of %i previously existing.\n',length(matchedIdxs),prevAvailableTracks);
                fprintf('Unmatched tracks: %i of %i previously existing.\n',length(unmatchedIdxs),prevAvailableTracks);
                fprintf('%i peaks remaining of %i.\n',size(peakMatrix,2),totalPeaks);
                fprintf('Deactivated tracks: %i\n',length(deactivatedIdxs)+length(deleteIdxs));
                fprintf('Deleted tracks: %i\n',length(deleteIdxs));
                activeIdxs = structArrayOperations(availableTracks,'status','==','active');
                asleepIdxs = structArrayOperations(availableTracks,'status','==','asleep');
                fprintf('Asleep tracks: %i\n',length(asleepIdxs));
                fprintf('Active tracks: %i\n',length(activeIdxs));
                fprintf('Existing tracks (asleep + active): %i\n',length(availableTracks));
                fprintf('Inactive tracks: %i\n',length(inactiveTracks));
                fprintf('Total tracks: %i\n',length(availableTracks)+length(inactiveTracks));
                if ((length(activeIdxs)+length(asleepIdxs))~=length(availableTracks))
                    error('Hmm... active + asleep is not equal to existing. Something is wrong and you dont know how to fix it...');
                end
                fprintf('Do all these make sense? asleep + active = existing is always necessary, so is existing + inactive = total.\n')
            end
        else
            if DEBUG == 1
                fprintf('There are no existing tracks or all tracks are inactive. Therefore no tracks were updated.\n');
                fprintf('This is frame: %i.\n',currentFrame);
                %fprintf('If this is not frame 1, then either your tracks are coming up too short or there is something wrong.\n\n');
            end
        end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Allocating new tracks -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-||

    % New tracks should be allowed if:  
        %   1) max length of a new track is greater than the minimum allowed track length.
        %   2) there is space for new tracks;
    if ((totalFrames-currentFrame+1) >= minTrackLength)
        % Calculating number of new tracks
        totalActiveTracks = length(structArrayOperations(availableTracks,'status','==','active')); % Gathering active indexes
        totalTracks = length(availableTracks)+length(inactiveTracks);
        totalNewTracks = min(size(peakMatrix,2),maxTracksPerFrame-totalActiveTracks); % Number of new tracks
        newTracks(1:totalNewTracks) = setNewTrack(); % Preallocating new tracks

        if DEBUG == 1
            fprintf('\nThis should allow for %i new tracks.\n',totalNewTracks);
        end

        newIdx = 1;
        while (newIdx <= totalNewTracks)
            newTracks(newIdx) = setNewTrack(peakMatrix(:,newIdx),currentFrame);
            totalTracks = totalTracks + 1;
            newIdx = newIdx + 1;
        end

        if DEBUG == 1
            fprintf('I have allocated %i tracks.\n',newIdx-1);
        end
        availableTracks = [availableTracks,newTracks];
    end

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Organizing track array -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
        currentTracks = [availableTracks,inactiveTracks];

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Treating last frame -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
        
        if currentFrame == totalFrames
            for trackIdx = 1:length(currentTracks)
                currentTracks(trackIdx) = setTrackInactive(currentTracks(trackIdx),currentFrame,totalFrames);
            end
        end

end

            