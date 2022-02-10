clear; clc; close all;

if isunix
    TFR_files_dir = '/Volumes/HD_MAURICIO/datasets/MedleyDB/MedleyDB_strc_tens_FChT_SLS_TFRs/';
else
    TFR_files_dir = '\Volumes\HD_MAURICIO\datasets\MedleyDB\MedleyDB_strc_tens_FChT_SLS_TFRs\';
end

file_name = files_in_path([TFR_files_dir 'FChT_ST_SLS'],{'mat'});

n_files = length(file_name);
% n_files = 200;

gini_all = [];

for ind_file = 1 : n_files

    clc;
    disp([num2str(100*ind_file/n_files) '%'])
%     disp(file_name(ind_file).clean_file_name)
    
    load([TFR_files_dir '/FChT_ST_SLS/' file_name(ind_file).clean_file_name]);
    TFR_FChT_ST_SLS = STFChT_comb;
    load([TFR_files_dir '/FChT_ST_LS/' file_name(ind_file).clean_file_name]);
    TFR_FChT_ST_LS = STFChT_comb;
    load([TFR_files_dir '/FChT_ST_SWGM/' file_name(ind_file).clean_file_name]);
    TFR_FChT_ST_SWGM = STFChT_comb;
    
    load([TFR_files_dir '/SLS/' file_name(ind_file).clean_file_name]);
    TFR_SLS = STFChT_comb;
    load([TFR_files_dir '/LS/' file_name(ind_file).clean_file_name]);
    TFR_LS = STFChT_comb;
    load([TFR_files_dir '/SWGM/' file_name(ind_file).clean_file_name]);
    TFR_SWGM = STFChT_comb;
    
    load([TFR_files_dir '/STFT-1024/' file_name(ind_file).clean_file_name]);
    load([TFR_files_dir '/STFT-2048/' file_name(ind_file).clean_file_name]);

    n_blocks = 10;
    n_frames_per_block = floor(size(TFR_SLS,2)/n_blocks);
    
    for block_ind = 1 : n_blocks
        gini = [];
        initial_frame = 1 + n_frames_per_block*(block_ind - 1);
        
        if block_ind < n_blocks
            final_frame = n_frames_per_block*block_ind;
        else
            final_frame = size(TFR_SLS,2);
        end

        gini(1, 1) = computeGiniIndex(TFR_FChT_ST_SLS(:, initial_frame:final_frame));
        gini(1, 2) = computeGiniIndex(TFR_FChT_ST_LS(:, initial_frame:final_frame));
        gini(1, 3) = computeGiniIndex(TFR_FChT_ST_SWGM(:, initial_frame:final_frame));
        gini(1, 4) = computeGiniIndex(TFR_SLS(:, initial_frame:final_frame));
        gini(1, 5) = computeGiniIndex(TFR_LS(:, initial_frame:final_frame));
        gini(1, 6) = computeGiniIndex(TFR_SWGM(:, initial_frame:final_frame));
        gini(1, 7) = computeGiniIndex(STFT_1024(:, initial_frame:final_frame));
        gini(1, 8) = computeGiniIndex(STFT_2048(:, initial_frame:final_frame));

        gini_all = cat(1, gini_all, gini);

    end
   

%     if ind_file == 1
%         [norm_samps_dist_1, ~] = compute_samps_dist(TFR);
%         [norm_samps_dist_2, ~] = compute_samps_dist(STFT_1024);
%         [norm_samps_dist_3, ~] = compute_samps_dist(STFT_2048);
%         [norm_samps_dist_4, ~] = compute_samps_dist(STFT_4096);
%         
%         samps_dist = [norm_samps_dist_1(1:samps_lims(end)), norm_samps_dist_2(1:samps_lims(end)),...
%                             norm_samps_dist_3(1:samps_lims(end)), norm_samps_dist_4(1:samps_lims(end))];
%         
%     elseif ind_file < n_files
%         [norm_samps_dist_1, ~] = compute_samps_dist(TFR);
%         [norm_samps_dist_2, ~] = compute_samps_dist(STFT_1024);
%         [norm_samps_dist_3, ~] = compute_samps_dist(STFT_2048);
%         [norm_samps_dist_4, ~] = compute_samps_dist(STFT_4096);
%         
%         samps_dist = cat(3, samps_dist, [norm_samps_dist_1(1:samps_lims(end)), norm_samps_dist_2(1:samps_lims(end)),...
%                                                       norm_samps_dist_3(1:samps_lims(end)), norm_samps_dist_4(1:samps_lims(end))]);
%     elseif ind_file == n_files
%         [norm_samps_dist_1, ~] = compute_samps_dist(TFR);
%         [norm_samps_dist_2, ~] = compute_samps_dist(STFT_1024);
%         [norm_samps_dist_3, ~] = compute_samps_dist(STFT_2048);
%         [norm_samps_dist_4, ~] = compute_samps_dist(STFT_4096);
%         
%         samps_dist = cat(3, samps_dist, [norm_samps_dist_1(1:samps_lims(end)), norm_samps_dist_2(1:samps_lims(end)),...
%                                                       norm_samps_dist_3(1:samps_lims(end)), norm_samps_dist_4(1:samps_lims(end)) ...
%                                                       ; zeros(size(samps_dist,1) - size(norm_samps_dist_1, 1) , 4)]);
% 
%     end
    
end

gini = gini_all;

% save('Statistics_DAFx', 'samps_dist');
save('Statistics_JAES', 'gini'); % Full spectrum
% save('Statistics_JAES_2', 'gini'); % Limited spectrum

%% 
clear; clc; close all;

load('Statistics_JAES.mat') % Full spectrum
% load('Statistics_JAES_2.mat') % Limited spectrum

if isunix
    figsPath = path_check('./figs/');
else
    figsPath = path_check('.\figs\');
end

% Rank with all TFRs
[~, Gini_rank] = sort(gini, 2, 'descend');
ranks = [sum(Gini_rank == 1) ;...
            sum(Gini_rank == 2) ;...
            sum(Gini_rank == 3) ;...
            sum(Gini_rank == 4) ;...
            sum(Gini_rank == 5) ;...
            sum(Gini_rank == 6) ;...
            sum(Gini_rank == 7) ;...
            sum(Gini_rank == 8)]*100/size(gini,1);
        
ranks = cumsum(ranks)

figure; b = bar(ranks)
legend('S-F-SLS', 'S-F-LS', 'S-F-SWGM', 'S-SLS', 'S-LS', 'S-SWGM', 'S(1024)', 'S(2048)',...
          'Location', 'NorthEastOutside')%, 'Orientation', 'horizontal')
xlabel('Rank Position')
ylabel('Top-X (%)')
axis([0.5 7.5 0 102])

tit = 'Statistics_JAES';

b(1).FaceColor = 'r';
b(2).FaceColor = 'b';
b(3).FaceColor = 'k';
b(4).FaceColor = 'c';
b(5).FaceColor = 'm';
b(6).FaceColor = 'g';
b(7).FaceColor = 'w';
b(8).FaceColor = 'y';

% b(1).LineStyle = '--';
% b(2).LineStyle = '--';
% b(3).LineStyle = '--';

%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 1200 900]); % Theses
%     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % Presentation
%     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 400 400]); % Small figures (3)
figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim', [1 1 800 300]);
figFileName = [figsPath tit];

formatFig(gcf, figFileName, 'en', figProp);


% % Boxplot normalized Gini
% gini_norm_factors = sum(gini,2);
% gini_norm_factors = repmat(gini_norm_factors,1,size(gini,2));
% 
% gini_norm = gini./gini_norm_factors;
% 
% TFR_names = {'SLS', 'LS', 'SWGM', 'S-1024', 'S-2048'};
% 
% figure; boxplot(gini, TFR_names, 'Notch','on')
% title('Gini Index')
% 
% figure; boxplot(gini_norm, TFR_names, 'Notch','on')
% title('Normalized Gini Index')

% % % Removing LS - - - - - - - - -
% gini_rm = gini;
% gini_rm(:, 2) = [];
% 
% % % Boxplot normalized Gini
% % gini_norm_factors = sum(gini_rm,2);
% % gini_norm_factors = repmat(gini_norm_factors,1,size(gini_rm,2));
% % 
% % gini_norm = gini_rm./gini_norm_factors;
% % 
% % TFR_names = {'SLS', 'SWGM', 'S-1024', 'S-2048'};
% % 
% % figure; boxplot(gini_norm, TFR_names, 'Notch','on')
% % title('Normalized Gini Index')
% % 
% % figure; hist(gini_norm)
% % title('Normalized Gini Index - No LS')
% 
% [~, Gini_rank] = sort(gini_rm, 2, 'descend');
% 
% % 1024 switched with 2048!!!
% figure; bar([sum(Gini_rank == 1) ; sum(Gini_rank == 2) ; sum(Gini_rank == 4); sum(Gini_rank == 3);] ...
%                                                                                                             *100/size(gini_rm,1))
% legend('SLS', 'SWGM', 'S-1024', 'S-2048', 'Location', 'NorthEastOutside')
% xlabel('Rank Positions')
% ylabel('Hits (%)')
% 
% tit = 'Statistics-JAES-no-LS';
% 
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 1200 900]); % Theses
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % Presentation
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 400 400]); % Small figures (3)
% figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim', [1 1 700 300]);
% figFileName = [figsPath tit];
% 
% formatFig(gcf, figFileName, 'en', figProp);

% 
% % Removing SWGM - - - - - - - - - - - -
% gini_rm = gini;
% gini_rm(:,3) = [];
% 
% [~, Gini_rank] = sort(gini_rm, 2, 'descend');
% 
% figure; bar([sum(Gini_rank == 1) ; sum(Gini_rank == 2) ; sum(Gini_rank == 3); sum(Gini_rank == 4);] ...
%                                                                                                             *100/size(gini_rm,1))
% legend('SLS', 'LS', 'S-1024', 'S-2048')
% xlabel('Rank Positions')
% ylabel('Hits (%)')
% 
% tit = 'Statistics-JAES-no-SWGM';
% 
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 1200 900]); % Theses
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % Presentation
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 400 400]); % Small figures (3)
% figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim', [1 1 800 400]);
% figFileName = [figsPath tit];
% 
% formatFig(gcf, figFileName, 'en', figProp);        
%         
% figure; bar(sum(Gini_rank == 1))
% title('Rank Positions')

% for samps_lim = samps_lims
%     
%     samps_lim
%     
%     samples_TFR = samps_dist(1:samps_lim, 1, :);
%     samples_S1 = samps_dist(1:samps_lim, 4, :);
%     samples_S2 = samps_dist(1:samps_lim, 3, :);
%     samples_S3 = samps_dist(1:samps_lim, 2, :);
% 
%     computeGiniIndex(samples_TFR(:))
%     computeGiniIndex(samples_S1(:))
%     computeGiniIndex(samples_S2(:))
%     computeGiniIndex(samples_S3(:))
%     disp('- - - - - -')
% 
%     TFR_names = {'C-TFR', 'S-4096', 'S-2048', 'S-1024'};
%     figure; boxplot([samples_TFR(:), samples_S1(:), samples_S2(:), samples_S3(:)], ...
%                 TFR_names, 'Notch','on')
% end

% tit = 'Statistics_DAFx';
% 
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
% %     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 1200 900]); % Theses
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % Presentation
% %     figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 400 400]); % Small figures (3)
% figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim', [1 1 700 500]);
% figFileName = [figsPath tit];
% 
% formatFig(gcf, figFileName, 'en', figProp);


%     for energy_ind = 1:length(energy_lims)
%         % Full band
%         [statistics(block_ind, 1, energy_ind, 1), ~] = compute_samps_percentage(uncompress(TFR_comb...
%                                     ((freq_combined<22500), initial_frame : final_frame)), energy_lims(energy_ind));        
%         [statistics(block_ind, 2, energy_ind, 1), ~] = compute_samps_percentage(uncompress(STFT_tensor...
%                                     ((freq_combined<22500), initial_frame : final_frame, 1)), energy_lims(energy_ind));
%         [statistics(block_ind, 3, energy_ind, 1), ~] = compute_samps_percentage(uncompress(STFT_tensor...
%                                     ((freq_combined<22500), initial_frame : final_frame, 2)), energy_lims(energy_ind));
%         [statistics(block_ind, 4, energy_ind, 1), ~] = compute_samps_percentage(uncompress(STFT_tensor...
%                                     ((freq_combined<22500), initial_frame : final_frame, 3)), energy_lims(energy_ind));        
%                                 
%         % Narrow band
%         [statistics(block_ind, 1, energy_ind, 2), ~] = compute_samps_percentage(uncompress(TFR_comb...
%                                     ((freq_combined<6000), initial_frame : final_frame)), energy_lims(energy_ind));        
%         [statistics(block_ind, 2, energy_ind, 2), ~] = compute_samps_percentage(uncompress(STFT_tensor...
%                                     ((freq_combined<6000), initial_frame : final_frame, 1)), energy_lims(energy_ind));
%         [statistics(block_ind, 3, energy_ind, 2), ~] = compute_samps_percentage(uncompress(STFT_tensor...
%                                     ((freq_combined<6000), initial_frame : final_frame, 2)), energy_lims(energy_ind));
%         [statistics(block_ind, 4, energy_ind, 2), ~] = compute_samps_percentage(uncompress(STFT_tensor...
%                                     ((freq_combined<6000), initial_frame : final_frame, 3)), energy_lims(energy_ind));
%                                 
%     end

%     % Full band
%     Gini_ind(block_ind, 1, 1) = computeGiniIndex(uncompress(TFR_comb...
%                                 ((freq_combined<22500), initial_frame : final_frame)));        
%     Gini_ind(block_ind, 2, 1) = computeGiniIndex(uncompress(STFT_tensor...
%                                 ((freq_combined<22500), initial_frame : final_frame, 1)));
%     Gini_ind(block_ind, 3, 1) = computeGiniIndex(uncompress(STFT_tensor...
%                                 ((freq_combined<22500), initial_frame : final_frame, 2)));
%     Gini_ind(block_ind, 4, 1) = computeGiniIndex(uncompress(STFT_tensor...
%                                 ((freq_combined<22500), initial_frame : final_frame, 3)));
% 
%     % Narrow band
%     Gini_ind(block_ind, 1, 2) = computeGiniIndex(uncompress(TFR_comb...
%                                 ((freq_combined<6000), initial_frame : final_frame)));        
%     Gini_ind(block_ind, 2, 2) = computeGiniIndex(uncompress(STFT_tensor...
%                                 ((freq_combined<6000), initial_frame : final_frame, 1)));
%     Gini_ind(block_ind, 3, 2) = computeGiniIndex(uncompress(STFT_tensor...
%                                 ((freq_combined<6000), initial_frame : final_frame, 2)));
%     Gini_ind(block_ind, 4, 2) = computeGiniIndex(uncompress(STFT_tensor...
%                                 ((freq_combined<6000), initial_frame : final_frame, 3)));