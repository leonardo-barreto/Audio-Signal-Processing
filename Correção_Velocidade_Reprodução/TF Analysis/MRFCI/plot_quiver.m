function ok = plot_quiver(varargin)
X = varargin{1};
Freq_orig = varargin{2}; 
Time_orig = varargin{3};
C = varargin{4};
thetas = varargin{5};

F_new = Freq_orig/20000;

time = Time_orig;

% figure;
% % subplot(1,2,1);
% imagesc(time,F,20*log10(S)); black = abs(1-gray); colormap(black);
% set(gca,'YDir','normal');
% hold on

% Vetores:
 [Taux,Faux] = meshgrid(time, F_new);

figure;
imagesc(time, F_new, X); black = abs(1-gray); colormap(black);
% imagesc(X); black = abs(1-gray); colormap(black);
set(gca,'YDir','normal');
hold on
eixoX = C.*cos(thetas);
eixoY = C.*sin(thetas);

quiver(Taux, Faux, eixoX, eixoY, .5, 'b', 'linewidth', 1.5); %axis equal
hold off

% axis([0.3 0.37 1.03 1.07]); % structure_tensor_zoom
% axis([0.24 0.40 1.01 1.09]); % curve
% axis([0.54 0.73 0.92 .98]); % curves
% axis([0.07 0.18 0.975 1.017]); % attack

% tit = 'ST_curves2';
% figProp = struct('size', .1,'font','Times','lineWidth',2,'figDim',[1 1 900 500]);
% figFileName = [figs_path tit];
% formatFig(gcf, figFileName, 'en', figProp);

ok = 1;

end