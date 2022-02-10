plot(xpdf,ypdf/sum(ypdf),'linewidth',1.5);
xlim([-4 4]);
xlabel('$\alpha$ values');
ylabel('Distribution');

figProp = struct('size', 32, 'font', 'Times', 'lineWidth', 2, 'figDim', [1 1 900 500]);
figFileName = [figs_path 'EstimatedHistogram'];
formatFig(gcf, figFileName, 'en', figProp);