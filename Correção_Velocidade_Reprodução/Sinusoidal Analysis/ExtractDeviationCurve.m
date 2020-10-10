function deviationCurve = ExtractDeviationCurve (frameArray,sinAnalysisParameters,signalTracks,DEBUG)

    fprintf('\n\n------- DEVIATION CURVE EXTRACTION STARTED ------\n\n');

    organizedTracks = sortStruct(signalTracks,'startFrame',1);
    totalTracks = length(signalTracks);

    totalFrames = frameArray(1).totalFrames;

    if DEBUG == 1
        fprintf('Signal has %i sinusoidal tracks.\n\n',length(organizedTracks));
    end

    trackFrequencies = {};
    deviationMatrix = NaN(totalTracks,totalFrames);
    trackMeanFrequencies = [];

    for trackIndex = 1:totalTracks
        trackFrequencies{trackIndex,1} = getfield(signalTracks(trackIndex),'frequencyEvolution',{[1:signalTracks(trackIndex).length]});
        trackMeanFrequencies(trackIndex) = mean(trackFrequencies{trackIndex});
    end

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

    %This calculates mean value ignoring NaNs, which 
    %are inexistent tracks. Therefore, it accounts only for existing tracks in a given frame.
        deviationCurve = [];
        deviationCurve(1:totalFrames) = nanmean(deviationMatrix(:,1:totalFrames)); 
    %Any NaN value means that there are no tracks in the frame, and therefore no deviation could be detected.
    %So, the resampling algorithm shouldn't change the original sampling rate (resampling factor of 1).
        deviationCurve(isnan(deviationCurve)) = 1; 

        

    % plot the deviation

    if DEBUG == 1  
        deviationPercentageMatrix = (deviationMatrix-1).*100; %This for plotting purposes only.  
        figure;
        surf(1:totalFrames,1:totalTracks,deviationPercentageMatrix);
        shading interp
        axis tight
        view(0, 90)
        xlabel('Frames')
        ylabel('Tracks (Indexes)')
        title('Deviation curve of each track around its own mean frequency');
        hcol = colorbar;
        ylabel(hcol, 'Deviation from mean frequency (%)');

        figure;
        plot(sinAnalysisParameters.timeInstants,deviationCurve,'LineWidth',3);
        X = sprintf('Curva de desvio relativo de pitch');
        title(X,'FontSize', 30);
        xlabel('Tempo (s)','FontSize', 30);
        ylabel('Frequencia relativa a media','FontSize', 30);
        set(gca,'FontSize', 30)
        %ylim([min(deviationCurve) max(deviationCurve)])
    end

    fprintf('------- DEVIATION CURVE EXTRACTED ------\n\n');

end