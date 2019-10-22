function organizedTracks = PlotTracks(frameArray,sinAnalysisParameters,trackArray);

    timeInstants = sinAnalysisParameters.timeInstants;

    organizedTracks = sortStruct(trackArray,'startFrame',1);

    figure;
    hold on;

    for trackIndex = 1:length(organizedTracks)

        trackStart = organizedTracks(trackIndex).startFrame;
        trackEnd = organizedTracks(trackIndex).finalFrame;
        trackFrames = trackStart:1:trackEnd;
        trackFrequencies = organizedTracks(trackIndex).frequencyEvolution;

        if (length(trackFrames) ~= length(trackFrequencies))
            fprintf('\nThis is track %i\n',trackIndex);
            error('Hmm. Wrong.');
        end

        timeRes = sinAnalysisParameters.timeInstants(10)-sinAnalysisParameters.timeInstants(9);
        trackTimes = 0:timeRes:sinAnalysisParameters.timeInstants(end)+2; %2 just a tolerance.
        trackTimes = trackTimes(trackStart:trackEnd);

        plot(trackTimes,trackFrequencies./1000,'LineWidth',2);

    end
    y = gca;
    set(y,'yscale','log')
    X = sprintf('Sinusoidal tracking');
    title(X);
    xlabel('Time(s)');
    ylabel('Frequency (kHz)');
    hold off;

end