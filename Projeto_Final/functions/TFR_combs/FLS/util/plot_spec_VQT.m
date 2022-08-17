function ok = plot_spec_VQT(tit, X, hop, fs, fmin, bins_per_octave, plot_par, figs_path)

black = abs(1-gray);

% spectr = spectr_orig(:,(t_orig>=plot_par(1) & t_orig<=plot_par(2)));
% t = t_orig(t_orig>=plot_par(1) & t_orig<=plot_par(2));

figure; imagesc((abs(flipud(X))));
xtickVec = 0:round(fs/hop):size(X,2)-1;
set(gca,'XTick',xtickVec);
ytickVec = 0:bins_per_octave:size(X,1)-1;
set(gca,'YTick',ytickVec);
ytickLabel = round(fmin * 2.^( (size(X,1)-ytickVec)/bins_per_octave));
set(gca,'YTickLabel',ytickLabel);
xtickLabel = 0 : length(xtickVec) ;
set(gca,'XTickLabel',xtickLabel);
xlabel('Time [s]', 'FontSize', 12, 'Interpreter','latex'); 
ylabel('Frequency [Hz]', 'FontSize', 12, 'Interpreter','latex');
set(gca, 'FontSize', 20);
colormap(black)

% grid on;
box on;
axis on;

figFileName = [figs_path tit];
xlabel('Time [s]'); ylabel('Frequency [Hz]');

figProp = struct('size', plot_par(5), 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
formatFig(gcf,figFileName,'en',figProp);

ok = 1;
end