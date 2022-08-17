function ok = plot_spec(tit, spectro, t, freq, plot_par, figs_path, redLines)

figure; imagesc(t, freq, spectro)
% title(tit);
colormap(1-gray)
axis(plot_par(1:4));
set(gca,'YDir','normal');
% colorbar

% To show onset and offset
if ~isempty(redLines)
    hold on;
    for ii = 1:length(redLines)
        y = 0:100:plot_par(4); 
        x = redLines(ii)*ones(size(y)); 
        plot(x,y,'--r');
    end
    hold off;
end

xlabel ('Time [s]');
ylabel ('Frequency [Hz]');
% pause;

figProp = struct('size',plot_par(5),'font','Times','lineWidth',2,'figDim',[1 1 900 500]);
figFileName = [figs_path tit];
formatFig(gcf,figFileName,'en',figProp);

save(figFileName, 'tit', 'spectro', 't', 'freq', 'plot_par');

ok = 1;
end

