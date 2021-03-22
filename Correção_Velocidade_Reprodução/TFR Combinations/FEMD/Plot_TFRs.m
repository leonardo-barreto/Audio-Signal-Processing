close all;

if isunix
    figsPath = path_check('./figs/');
else
    figsPath = path_check('.\figs\');
end

axis_vec = [0 1 0 4];


% % LS - - - - - - - - - - - - - - - - - - - - 
% load(['./TFRs/FChT_ST_LS/' signal_name])

 
% figure; imagesc(time_2, F_final/1000, compress(STFChT_comb));
% set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
% % title('LS')
% % xlabel('Time (s)')
% % ylabel('Frequency (kHz)')
% axis(axis_vec)
% set(gca,'YTick',[]);
% set(gca,'XTick',[]);
% % pos = get(gca, 'Position');
% % pos(1) = pos(1) - .1; 
% % pos(3) = pos(3) - .1;
% % set(gca, 'Position', pos)
% 
% tit = ['LS_zoom'];
% 
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 1200 900]); % Theses
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % Presentation
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 400 400]); % Small figures (3)
% figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim', [1 1 400 400]);
% figFileName = [figsPath tit];
% 
% formatFig(gcf, figFileName, 'en', figProp);
% 

% SLS - - - - - - - - - - - - - - - - - - - - 
load(['./TFRs/FChT_ST_SLS/' signal_name])

figure; imagesc(time_2, F_final/1000, compress(STFChT_comb));
set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
% title('SLS')
xlabel('Time (s)')
ylabel('Frequency (kHz)')
axis(axis_vec)
% set(gca,'YTick',[]);
% set(gca,'XTick',[]);
% pos = get(gca, 'Position');
% pos(1) = pos(1) - .05; 
% pos(3) = pos(3) - .1;
% set(gca, 'Position', pos)

tit = 'SLS_no_filter';

%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 1200 900]); % Theses
%     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % Presentation
%     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 400 400]); % Small figures (3)
figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim', [1 1 800 300]);
figFileName = [figsPath tit];

formatFig(gcf, figFileName, 'en', figProp);

% 
% 
% % SWGM - - - - - - - - - - - - - - - - - - - - 
% load(['./TFRs/FChT_ST_SWGM/' signal_name])
% 
% % figure; 
% subplot(1,3,3); imagesc(time_2, F_final/1000, compress(STFChT_comb));
% set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
% title('SWGM')
% xlabel('Time (s)')
% % ylabel('Frequency (kHz)')
% axis(axis_vec)
% set(gca,'YTick',[]);
% % set(gca,'XTick',[]);
% pos = get(gca, 'Position');
% pos(1) = pos(1) - .1; 
% % pos(3) = pos(3) - .2;
% set(gca, 'Position', pos)




% % STFT-1024
% load(['./TFRs/STFT-1024/' signal_name])
% 
% figure; imagesc(time_2, F_final/1000, compress(STFT_1024));
% set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
% title('STFT-1024')
% xlabel('Time (s)')
% ylabel('Frequency (kHz)')
% axis(axis_vec)
% 
% % STFT-2048
% load(['./TFRs/STFT-2048/' signal_name])
% 
% figure; imagesc(time_2, F_final/1000, compress(STFT_2048));
% set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
% title('STFT-2048')
% xlabel('Time (s)')
% ylabel('Frequency (kHz)')
% axis(axis_vec)
