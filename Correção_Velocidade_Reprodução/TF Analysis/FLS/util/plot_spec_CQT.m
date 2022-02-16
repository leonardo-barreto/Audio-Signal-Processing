function ok = plot_spec_CQT(tit, spectr_orig, t_orig, freq_bins_CQT, Xcqt, plot_par, figs_path)

spectr = spectr_orig(:,(t_orig>=plot_par(1) & t_orig<=plot_par(2)));
t = t_orig(t_orig>=plot_par(1) & t_orig<=plot_par(2));

figure; surf(t,1:size(spectr,1),spectr,'EdgeColor','none');

bins = Xcqt.bins;
octaveNr = Xcqt.octaveNr;
fmin = Xcqt.fmin;

figFileName = [figs_path tit];
% save(figFileName,'tit', 'spectr_orig', 't_orig', 'freq_bins_CQT', 'bins', 'octaveNr', 'fmin');
caxis([min(min(spectr)) max(max(spectr))])
set(gca,'YTick',1:bins:octaveNr*bins);
h = get(gca); yTick = h.YTick';
yTickLabel = num2str(round(fmin*2.^((yTick-1)/bins)),5);

axis('tight'); view(0,90);
colormap(1-gray)
% colormap(gray)

set(gca,'YTickLabel',yTickLabel);
xlabel('Time [s]'); ylabel('Frequency [Hz]');

grid on;
box on;
axis on;

figProp = struct('size',plot_par(5),'font','Times','lineWidth',2,'figDim',[1 1 900 500]);
formatFig(gcf,figFileName,'en',figProp);
% close;

ok = 1;
end