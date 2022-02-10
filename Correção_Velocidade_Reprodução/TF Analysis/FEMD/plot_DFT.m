% Plots the spectrogram saved

% plot_par = [t_i, t_f, f_i, f_f]
plot_par = [3, 8, 0, 8000];

black = abs(1-gray);
imagesc(t, freq, spectrogram)

% Formating and saving the figure
title(tit);
colormap(black)
axis(plot_par(1:4));
set(gca,'YDir','normal');

xlabel ('Time [s]');
ylabel ('Frequency [Hz]');
title(tit)

figProp = struct('size', 35,'font','Times','lineWidth',2,'figDim',[1 1 900 500]);
figFileName = [figs_path tit];
formatFig(gcf, tit, 'en', figProp);