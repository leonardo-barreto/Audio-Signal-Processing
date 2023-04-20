function organizedTracks = PlotTracks(trackArray,timeInstants);

    if (~size(trackArray,2))
        organizedTracks = [];
        fprintf('\nWell... no tracks.\n');
        return;
    end

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

       
    plot(trackTimes,trackFrequencies,'LineWidth',4,'Color','k');
    plot(0,0,'Color','w');

    end
    y = gca;
    y_lim_vec = [0 5000];
    x_lim_vec = [0 timeInstants(end)];
    %set(y,'yscale','log')
    set(y,'YDir','normal');    
    ylim(y_lim_vec); 
    xlabel('Tempo (s)','FontSize', 30);
    ylabel('Frequencia (Hz)','FontSize', 30);
    %set(y,'FontSize', 35)
    %figProp = struct('size', 35, 'font', 'Times', 'lineWidth', 4 , 'figDim', [1 1 900 500]);
    %figFileName = 'Tracking';
    %formatFig(gcf, figFileName, 'en', figProp);
    hold off;

end