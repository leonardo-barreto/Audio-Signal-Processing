function [test1,test2] = ExtractDeviationCurve (frameArray,signalTracks)

    organizedTracks = sortStruct(signalTracks,'startFrame',1);
    totalTracks = length(signalTracks);

    totalFrames = frameArray(1).totalFrames;

    fprintf('\nFound %i tracks.\n',length(organizedTracks));

    trackFrequencies = {};
    deviationMatrix = zeros(totalTracks,totalFrames);
    trackMeanFrequencies = [];

    for trackIndex = 1:totalTracks
        trackFrequencies{trackIndex,1} = getfield(signalTracks(trackIndex),'frequencyEvolution',{[1:signalTracks(trackIndex).length]});
    end

    for trackIndex = 1:totalTracks
        trackMeanFrequencies(trackIndex) = mean(trackFrequencies{trackIndex});
    end

    test1 = trackMeanFrequencies;

    for frameIndex = 1:totalFrames

        tracksAbove = structArrayOperations(signalTracks,'startFrame','<=',frameIndex);
        tracksBelow = structArrayOperations(signalTracks,'finalFrame','>=',frameIndex);
        activeTracks = intersect(tracksAbove,tracksBelow);

        trackFrequency = [];
        for trackIndex = activeTracks
            trackStart = signalTracks(trackIndex).startFrame;
            trackFrequency = getfield(signalTracks(trackIndex),'frequencyEvolution',{[frameIndex+1-trackStart]});
            deviationMatrix(trackIndex,frameIndex) = trackFrequency/trackMeanFrequencies(trackIndex);
        end



    
    end

    test2 = deviationMatrix;

end