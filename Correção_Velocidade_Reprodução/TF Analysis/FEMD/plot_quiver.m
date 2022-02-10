F_new=F/20000;

% figure;
% % subplot(1,2,1);
% imagesc(time,F,20*log10(S)); black = abs(1-gray); colormap(black);
% set(gca,'YDir','normal');
% hold on

% Vetores:
 [Taux,Faux] = meshgrid(time,F_new);

figure;
imagesc(time, F_new, X.^.25); black = abs(1-gray); colormap(black);
set(gca,'YDir','normal');
hold on
eixoX = C.*cos(theta);
eixoY = C.*sin(theta);

quiver(Taux, Faux, eixoX, eixoY, .5, 'r', 'linewidth', 1); %axis equal
hold off
% axis([.07 .15 .73 .79])
tit = 'theta_arrows';

% figProp = struct('size', .1,'font','Times','lineWidth',2,'figDim',[1 1 400 300]);
% figFileName = [figs_path tit];
% formatFig(gcf, figFileName, 'en', figProp);


% figure;
% imagesc(time, F, X.^.25); black = abs(1-gray); colormap(black);
% set(gca,'YDir','normal');
% hold on
% eixoX = cos(alphas);
% eixoY = sin(alphas);
% quiver(Taux, Faux, eixoX, eixoY, .3, 'r', 'linewidth', 2); %axis equal
% hold off
% % axis([.23 .46 1.4 2])
% title('Usando somente o \alpha')

% figure;
% imagesc(time,F,C);
% colormap(black);
% set(gca,'YDir','normal');
% title('C')

% figure;
% plot(mean(alphas))
