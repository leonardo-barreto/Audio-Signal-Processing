function currentTracks = PartialTracking_2023MQ(inputFrame,totalFrames,currentTracks,DEBUG)

    DEBUG = 0;

    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Defining parameters -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|

        maxTracksPerFrame = 100;

        freqTolerance = (power(2,1/24)-1); % about 3% (quarter-tone)
        powerTolerance = 3;                % in dB
        maxHysteresis = 3;                 % in frames
        minTrackLength = 10;               % in frames
        maxTrackFrequency = 5000;          % in Hz
        minTrackPower = -60;               % in dB

        %POWER in ROW 1, FREQUENCY in ROW 2
            powerRow = 1;
            freqRow = 2;
        % Gathering frame data
            currentFrame = inputFrame.currentFrame;
            peakMatrix = inputFrame.peakMatrix;
            lastFrame = totalFrames;

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
        % Creating necessary arrays
        matchedIdxs = [];
        unmatchedIdxs = [];
        deactivatedIdxs = [];
        % Ordering tracks in frequency
        availableTracks = sortStruct(currentTracks, 'currentFrequency');
        if DEBUG == 1
            fprintf('Available (active or asleep) tracks: %i of %i.\n',length(availableTracks),length(currentTracks));
            fprintf('Existing peaks: %i\n',totalPeaks);
            fprintf('Proceeding to match tracks to peaks...\n\n');
        end 
        
    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Existing Tracks processing -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
        
        % Won't attempt to process existing tracks if there aren't any.
        if (~isempty(availableTracks))
            for trackIdx = 1:length(availableTracks)
                trackFrequency = availableTracks(trackIdx).frequencyEvolution(end);
                matchedPeakIdxs = intersect(find(peakMatrix(freqRow,:) < trackFrequency*(1 + freqTolerance)),find(peakMatrix(freqRow,:) > trackFrequency*(1 - freqTolerance)));
                %Power matching criterion (experimental)
                %matchedPeaksPower = intersect(find(peakMatrix(powerRow,:) < availableTracks(trackIdx).powerEvolution(end) + powerTolerance),find(peakMatrix(powerRow,:) > availableTracks(trackIdx).powerEvolution(end) - powerTolerance));
                %matchedPeakIdxs = intersect(matchedPeakIdxs,matchedPeaksPower);

                if (~isempty(matchedPeakIdxs)) % 1) Track matched to a peak
                    [freqGap,peakIdx] = min(abs(peakMatrix(freqRow,matchedPeakIdxs)-trackFrequency));
                    peakMatchIdx = matchedPeakIdxs(peakIdx);
                    if DEBUG == 1
                        fprintf('Peak %i matched to track %i.\n',peakMatchIdx,trackIdx);
                    end
                    % 2) Check if there is better match to the peak
                    betterMatch = find(abs(peakMatrix(freqRow,peakMatchIdx)-getArrayFields(availableTracks(trackIdx+1:end),'currentFrequency')) < freqGap)
                    if (isempty(betterMatch)) % current track is best match
                        availableTracks(trackIdx) = setTrackActive(availableTracks(trackIdx),peakMatrix(:,peakMatchIdx),currentFrame);
                        matchedIdxs(end+1) = trackIdx;
                        peakMatrix(:,peakMatchIdx) = [];
                        if DEBUG == 1
                            fprintf('Peak %i definite match to track %i.\n',peakMatchIdx,trackIdx);
                        end  
                    else % 2.1) Better match for peak found
                        if DEBUG == 1
                            fprintf('Better match for peak %i found.\n',peakMatchIdx,betterTrackIdx);
                        end
                        if (numel(matchedPeakIdxs) > 1 && peakMatchIdx > 1) % 2.1.a) Track can still be matched to adjacent peak
                            availableTracks(trackIdx) = setTrackActive(availableTracks(trackIdx),peakMatrix(:,peakMatchIdx-1),currentFrame);
                            matchedIdxs(end+1) = trackIdx;
                            peakMatrix(:,peakMatchIdx-1) = [];
                            if DEBUG == 1
                                fprintf('Track %i matched definitely to peak %i instead.\n',trackIdx,peakMatchIdx-1);
                            end
                        else % 2.1.b) Track can't be matched
                            unmatchedIdxs(end+1) = trackIdx;
                            if DEBUG == 1
                                fprintf('No definite match possible for track %i.\n',trackIdx);
                            end
                        end 
                    end 
                else % 1) No matches found
                    unmatchedIdxs(end+1) = trackIdx;
                    if DEBUG == 1
                        fprintf('Track %i NOT matched to any peak.\n',trackIdx);
                    end
                end 
            end
        end
        for trackIdx = unmatchedIdxs
            switch availableTracks(trackIdx).status
                case 1 %Track is active and must go to sleep 
                    availableTracks(trackIdx) = setTrackAsleep(availableTracks(trackIdx),currentFrame);
                case 2 %Track is asleep
                    if (availableTracks(trackIdx).hysteresis >= maxHysteresis) %Track is over hysteresis limit and must be deactivated.
                        availableTracks(trackIdx) = setTrackInactive(availableTracks(trackIdx),currentFrame,lastFrame);
                        deactivatedIdxs(end+1) = trackIdx;
                    else %Track is asleep and can continue to be so.
                        availableTracks(trackIdx) = setTrackAsleep(availableTracks(trackIdx),currentFrame);
                    end
            end
        end

        % TO DO: 
        %       CHECK UNMATCHED IDXS TREATMENT
        %       DELETE SHORT TRACKS
        %       NEW TRACKS
        %       ORGANIZE TRACK ARRAY FOR FUNCTION RETURN
        %       TREAT LAST FRAME (DEACTIVATE ALL TRACKS)

            