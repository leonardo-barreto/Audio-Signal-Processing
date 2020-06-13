function organizedTracks = PlotTracks(frameArray,sinAnalysisParameters,trackArray);

    if (~size(trackArray,2))
        organizedTracks = [];
        fprintf('\nWell... no tracks.\n');
        return;
    end

    timeInstants = sinAnalysisParameters.timeInstants;

    organizedTracks = sortStruct(trackArray,'startFrame',1);

    figure;
    hold on;

    colorVector = ['r-','b-','g-','m-','k-'];


    for trackIndex = 1:length(organizedTracks)

        trackStart = organizedTracks(trackIndex).startFrame;
        trackEnd = organizedTracks(trackIndex).finalFrame;
        trackFrequencies = organizedTracks(trackIndex).frequencyEvolution;
        trackTimes = timeInstants(trackStart:trackEnd);

        if (length(trackTimes) ~= length(trackFrequencies))
            fprintf('\nThis is track %i\n',trackIndex);
            fprintf('%i Track Times and %i trackFrequencies',length(trackTimes), length(trackFrequencies));
            error('Hmm. Wrong.');
        end

        

    plot(trackTimes,trackFrequencies,'LineWidth',5);

    end
    y = gca;
    set(y,'yscale','log')
    X = sprintf('Rastreamento de trilhas senoidais');
    title(X);
    xlabel('Tempo (s)','FontSize', 30);
    ylabel('Frequencia (Hz)','FontSize', 30);
    set(y,'FontSize', 30)
    %ylim([1997 2001])
    hold off;

end