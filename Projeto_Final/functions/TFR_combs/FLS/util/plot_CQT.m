% Plots the spectrogram saved

% plot_par = [t_i, t_f, f_i, f_f]
plot_par = [3, 8, 0, 8000];

black = abs(1-gray);

spectr = spectr_orig(:,(t_orig>=plot_par(1) & t_orig<=plot_par(2)));
t = t_orig(t_orig>=plot_par(1) & t_orig<=plot_par(2));

surf(t,1:size(spectr,1),spectr,'EdgeColor','none');

caxis([min(min(spectr)) max(max(spectr))])
axis('tight'); view(0,90);
set(gca,'YTick',1:bins:octaveNr*bins);
h = get(gca); yTick = h.YTick';
yTickLabel = num2str(round(fmin*2.^((yTick-1)/bins)),5);
colormap(black)
% title(tit)

set(gca,'YTickLabel',yTickLabel);
xlabel('Time [s]'); ylabel('Frequency [Hz]');

box on;
axis on;

figProp = struct('size', 35, 'font','Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
formatFig(gcf, tit, 'en',figProp);