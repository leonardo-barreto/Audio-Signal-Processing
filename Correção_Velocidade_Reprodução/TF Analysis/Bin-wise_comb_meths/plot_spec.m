function ok = plot_spec(tit, spectro, t, freq, plot_par, figs_path)

figure;
black = abs(1-gray);
imagesc(t, freq, spectro)

colormap(black)
axis(plot_par(1:4));
set(gca,'YDir','normal');

xlabel ('Time [s]');
ylabel ('Frequency [Hz]');

figProp = struct('size',plot_par(5),'font','Times','lineWidth',2,'figDim',[1 1 900 500]);
figFileName = [figs_path tit];
formatFig(gcf,figFileName,'en',figProp);

save(figFileName, 'tit', 'spectro', 't', 'freq', 'plot_par');

ok = 1;
end

