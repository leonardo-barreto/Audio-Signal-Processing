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
        % Ordering tracks in frequency
        availableTracks = sortStruct(currentTracks, 'currentFrequency');
        if DEBUG == 1
            fprintf('Available (active or asleep) tracks: %i of %i.\n',length(availableTracks),length(currentTracks));
            fprintf('Existing peaks: %i\n',totalPeaks);
            fprintf('Proceeding to match tracks to peaks...\n\n');
        end 
        
    % -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-| Matching tracks to peaks -|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|
        
        % Won't attempt to process existing tracks if there aren't any.
        if (~isempty(availableTracks))
            for trackIdx = 1:length(availableTracks)
                trackFrequency = availableTracks(trackIdx).frequencyEvolution(end);
                matchedPeaks = intersect(find(peakMatrix(freqRow,:) < trackFrequency*(1 + freqTolerance)),find(peakMatrix(freqRow,:) > trackFrequency*(1 - freqTolerance)));
                %Power matching criterion (experimental)
                %matchedPeaksPower = intersect(find(peakMatrix(powerRow,:) < availableTracks(trackIdx).powerEvolution(end) + powerTolerance),find(peakMatrix(powerRow,:) > availableTracks(trackIdx).powerEvolution(end) - powerTolerance));
                %matchedPeaks = intersect(matchedPeaks,matchedPeaksPower);

                if (~isempty(matchedPeaks)) % Track matched to a peak
                    [freqGap,matchIdx] = min(abs(peakMatrix(freqRow,matchedPeaks)-trackFrequency));
                    matchIdx = matchedPeaks(matchIdx);
                    trackFrequencies = getArrayFields(availableTracks,'currentFrequency');
                    trackFrequencies = trackFrequencies(1:end ~= matchIdx);
                    matchTracks = intersect(find(availableTracks),find());

                    % Old criterion
                    %availableTracks(trackIdx) = setTrackActive(availableTracks(trackIdx),peakMatrix(:,matchIdx),currentFrame);
                    %matchedIndexes(end+1) = trackIdx;
                    %peakMatrix(:,matchIdx) = []; % Removes from selection a peak that was used
                    if DEBUG == 1
                        fprintf('Peak %i matched to track %i.\n',matchIdx,trackIdx);
                    end
                else % No matches found
                    unmatchedIndexes(end+1) = trackIdx;
                    switch availableTracks(trackIdx).status
                        case 1 %Track is active and must go to sleep 
                            availableTracks(trackIdx) = setTrackAsleep(availableTracks(trackIdx),currentFrame);
                        case 2 %Track is asleep
                            if (availableTracks(trackIdx).hysteresis >= maxHysteresis) %Track is over hysteresis limit and must be deactivated.
                                availableTracks(trackIdx) = setTrackInactive(availableTracks(trackIdx),currentFrame,lastFrame);
                                deactivatedIndexes(end+1) = trackIdx;
                            else %Track is asleep and can continue to be so.
                                availableTracks(trackIdx) = setTrackAsleep(availableTracks(trackIdx),currentFrame);
                            end
                    end
                    if DEBUG == 1
                        fprintf('Track %i NOT matched to any peak.\n',trackIdx);
                    end
                end 
            end
        end    

            