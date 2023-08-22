function variationCurve = ExtractPVC(signalTracks,totalFrames)

    DEBUG = 0;

    % Gathering tracks and info
    organizedTracks = sortStruct(signalTracks,'startFrame',1);
    totalTracks = length(signalTracks);

    if DEBUG == 1
        fprintf('Signal has %i sinusoidal tracks.\n\n',length(organizedTracks));
    end

    % Calculating track mean frequencies
    trackFreqs = cell(totalTracks,1);
    trackMeanFreqs = zeros(totalTracks,1);
    for trackIdx = 1:totalTracks
        trackFreqs{trackIdx} = getfield(signalTracks(trackIdx),'frequencyEvolution',{[1:length(signalTracks(trackIdx).frequencyEvolution)]});
        trackMeanFreqs(trackIdx) = mean(trackFreqs{trackIdx});
    end

    % Calculating percentual frequency deviation on each track, per frame
    deviationMatrix = NaN(totalTracks,totalFrames);
    for frameIdx = 1:totalFrames
        tracksAbove = structArrayOperations(signalTracks,'startFrame','<=',frameIdx);
        tracksBelow = structArrayOperations(signalTracks,'finalFrame','>=',frameIdx);
        activeTracks = intersect(tracksAbove,tracksBelow);
        trackFreq = [];
        for trackIdx = activeTracks
            trackFreq = getfield(signalTracks(trackIdx),'frequencyEvolution',{[frameIdx+1-signalTracks(trackIdx).startFrame]});
            deviationMatrix(trackIdx,frameIdx) = trackFreq/trackMeanFreqs(trackIdx);
        end
    end

    %This calculates mean value ignoring NaNs, which 
    %are inexistent tracks. Therefore, it accounts only for existing tracks in a given frame.
        variationCurve = [];
        variationCurve(1:totalFrames) = nanmean(deviationMatrix(:,1:totalFrames)); 
    %Any NaN value means that there are no tracks in the frame, and therefore no deviation could be detected.
    %So, the resampling algorithm shouldn't change the original sampling rate (resampling factor of 1).
        variationCurve(isnan(variationCurve)) = 1;
end