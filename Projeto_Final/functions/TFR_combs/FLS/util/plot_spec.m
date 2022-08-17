function ok = plot_spec(tit, spectro, t, freq, plot_par, figs_path, redLines)
%PLOT_SPEC Summary of this function goes here
%   Detailed explanation goes here

figure;
black = abs(1-gray);
imagesc(t, freq, spectro)
% title(tit);
colormap(black)
axis(plot_par(1:4));
% axis([0, .8, 0, 2000])
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
% close;

% save(figFileName, 'tit', 'spectro', 't', 'freq', 'plot_par');

ok = 1;
end

