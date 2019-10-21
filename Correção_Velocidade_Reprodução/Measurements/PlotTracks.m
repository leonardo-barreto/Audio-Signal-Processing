function organizedTracks = PlotTracks(trackArray);

    organizedTracks = sortStruct(trackArray,'startFrame',1);

    fprintf('\nFound %i tracks.\n',length(organizedTracks));

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

        plot(trackFrames,trackFrequencies,'LineWidth',2);

    end
    y = gca;
    set(y,'yscale','log')
    X = sprintf('Sinusoidal tracking');
    title(X);
    xlabel('Frames');
    ylabel('Frequency (Hz)');
    hold off;

end