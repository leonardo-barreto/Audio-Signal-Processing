%clear variables; clc; close all;

%%%%%%%%%%%%%
%       Setting Paths  
%%%%%%%%%%%%%

figs_path = path_check('./figs_DSc/'); 

if isunix
    addpath('./../../Audio/sinais_selecionados')
    addpath('./../../Audio')
    addpath('./Sinais')
    addpath('./read_signals')
%     addpath('./Resultados Saliencia')
else
    addpath('.\..\..\Audio')
    addpath('.\..\..\Audio\sinais_selecionados')
    addpath('.\Sinais')
    addpath('.\read_signals')
%     addpath('./Resultados Saliencia')
end

%%%%%%%%%%%%%%%%%%
%       Parameters configuration     %
%%%%%%%%%%%%%%%%%%

% Analysis TFR - - - - - - - - - - - - - - - - - - - -
fs = 44100;
Nf = 1024; hop = Nf/4;
comb_method = 3; % 1- SWGM, 2 - LS, 3 - SLS, 4 - All methods
combine_original_TFR = 1; % Includes the initial TFR used by the ST for the alphas estimation <<<<<<<<<<<<<
combine_large_STFT = 1; % Includes the STFT computed with the same window size as the STFChT
combine_extra_large_STFT = 0; % Includes an STFT computed with a larger window size than the STFChT
Nf_3 = Nf*8;

% use_init_comb = 0; % For combining FChT instances before the ST procedure
% initial_alphas = [-3 3]; % For combining FChT instances before the ST procedure

% Final TFR - - - - - - - - - - - - - - - - - - -
resMult = 2; % Multiplies the resolution of the Analysis TFR
Nf_2 = Nf*resMult;
hop_2 = hop;

% useTheoAlphas = 0;  % For testing with controlled signals
% plotTheoAlphas = 0; % For testing with controlled signals
renyi_ord = 3;% For evaluating the sparsity of final results
Num_it = 1;   % Number of iterations
filter_alphas_flag = 1; % Alphas filtering
ord = 5;% Alphas filtering
ord_median = 5;% Alphas filtering

% - - - Structure Tensor - - -
sigma_n = 1.5; sigma_k = 1.5; Ngauss = 9; % G parameters
h = 0.1; % Standard deviation for the distribution
percent_trsh = 0.05; % Percent of the highest peak to consider as relevant (GMM)
epsilon = 1;

% - - - Local Sparsity Combination - - -
switch comb_method
    case 2 % LS
        zeta = 202; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
        eta = 1; % controls the energy weighting process
    case 3 % SLS
        zeta = 70; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
        eta = 1; % controls the energy weighting process
end

time_span = 0.03; % time span of the 2D window in seconds
Ncols_DFT = ceil(fs*time_span/hop);
size_w_S = [21 Ncols_DFT];
size_w_E = [9  2*Ncols_DFT];

% plot_window_E(size_w_E, 'window_E', 27)
% plot_window_S(size_w_S, 'window_S', 27)

% return

% - - - SWGM Combination - - -
beta = .5;

% - - - Plotting parameters - - -
plot_par = [0.05, .95, 0, 8000, 35];
ylim_vec = [12 19]; % in kHz
ylim_alphas = [];

%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%%%%%%%%%%%%%%%%
%       Reading input signal      %
%%%%%%%%%%%%%%%%

% read_harmonic_sin_1;
% read_violins_vibrato;

% read_violins_vibrato_drums; % < < < < < < < < < < <
% read_orchestra_vocal_drum % < < < < < < < < < < <
% read_pop; % < < < < < < < < < < <
% read_opera_vocals;
% read_operaTenor;
% read_violins_piano; % xxxxx
% read_signals; % xxxxx

%%% For AES Journal %%%
% read_harmonic_sin; % < < < < < < < < < < <
% read_violins_drum2; % < < < < < < < < < < <
% read_violin_vocal; % < < < < < < < < < < <
% read_violins_vibrato

%%% For Thesis %%%
read_harmonic_sin; % < < < < < < < < < < <
% read_violin_vocal; % < < < < < < < < < < <
% read_violins_drum2; % < < < < < < < < < < <

% TFR_files_dir = './TFRs/';

% % - - -  - - - - - - - - - - - - - - - - - - - - MedleyDB - - - - - - - - - - - - - - - - - - - - - - - - - -
% %Get file names (wav)
% audio_files_dir = '/Volumes/HD_MAURICIO/datasets/MedleyDB/MedleyDB_V1_mixes_cut_10';
% wav_file_name = files_in_path(audio_files_dir,{'wav'});
% 
% TFR_files_dir = '/Volumes/HD_MAURICIO/datasets/MedleyDB/MedleyDB_strc_tens_FChT_SLS_TFRs';
% 
% NumOfSources = 2;
% useTheoAlphas = 0; % For testing with controlled signals
% plotTheoAlphas = 0; % For testing with controlled signals
% 
% % for file_ind = 35:length(wav_file_name)
% file_ind = 1;
% 
% %     tic;
%     
%     %Load wav file - - - - - - - - - - - - - - - - - - - - -
%     [data, fs_orig] = audioread(wav_file_name(file_ind).name);
%     
%     signal_name = wav_file_name(file_ind).clean_file_name
%         
%     % Make it mono
%     if size(data,2) > 1
%         x = mean(data.');
%     else
%         x = data;
%     end
%     x = x(:);
%     
%     % Resampling
%     if fs_orig ~= fs
%         x = resample(x,fs,fs_orig);
%     end    
%     % - - - - - - - - - - - - - - - - - - - - - - - MedleyDB - - - - - - - - - - - - - - - - - - - - - - - - - -
%     
    % soundsc(x,fs)
    % Normalizing x
    x = x/sqrt(mean(x.^2));

    %% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    %%%%%%%%%%%%%%
    %       Structure Tensor     %
    %%%%%%%%%%%%%%

    % Computing the STFT:------------------------------------------------------
    tic_all = tic;
    [STFT,F,time] = spectrogram(x, hanning(Nf),Nf-hop,Nf,fs); 

    X_original = abs(STFT); 

    % % Using initial combination of STFChTs:---------------------------------------
    % if use_init_comb
    %     Num_frames = length(time); STFChT_init = zeros(Nf/2+1, Num_frames, length(initial_alphas));
    %     for ii = 1:Num_frames
    %         for alphas_ind = 1:length(initial_alphas)
    %             t0 = (ii-1)*hop+1;
    %             xn = x(t0:t0+Nf-1);
    %             [~,mag_specs] = calculaFChT(xn, fs, Nf, 1, initial_alphas(alphas_ind));
    %             STFChT_init(:,ii,alphas_ind) = abs(mag_specs(1:Nf/2+1));
    %         end
    %     end
    %     % Combining the representations
    %     initial_spectrograms = cat(3, STFChT_init, X_original);
    %     X_original = spectrogram_comb_local_sparsity(initial_spectrograms, size_w_S, size_w_E, zeta, eta);
    %     X = X_original;
    % else
        X = X_original;
    % end

    % % Plotting - - - -
    % figure;
    % imagesc(time,F,20*log10(S))
    % colormap(black); set(gca,'YDir','normal');
    % ylim(ylim_vec)
    % title('Initial combination')
    % % Plotting - - - -

    % Computing the structure tensor:------------------------------------------

    % Using imgradientxy:
    [Xn, Xw] = imgradientxy(compress(X));

    %     %Plotting - - - -
    %     figure; subplot(1,2,1); imagesc(time,F,Sb); 
    %     colormap(black); set(gca,'YDir','normal'); ylim([0 20000])
    %     title('Time Derivative','FontName','Times','FontSize',14)
    %     subplot(1,2,2); imagesc(time,F,Sk); 
    %     colormap(black); set(gca,'YDir','normal'); ylim([0 20000])
    %     title('Frequency Derivative','FontName','Times','FontSize',14)
    %     %Plotting - - - -

    % Computing T11, T21, T12 and 22:
    G = gauss2D(Ngauss, [sigma_n sigma_k]);
    T11 = conv2(Xn.*Xn, G, 'same');
    T21 = conv2(Xw.*Xn, G, 'same'); T12 = T21;
    T22 = conv2(Xw.*Xw, G, 'same');

    %     %Plotting - - - -
    %     figure; subplot(1,3,1); imagesc(time,F,T11); 
    %     colormap(black); set(gca,'YDir','normal'); ylim([0 20000])
    %     title('T11','FontName','Times','FontSize',14)
    %     subplot(1,3,2); imagesc(time,F,T21); 
    %     colormap(black); set(gca,'YDir','normal'); ylim([0 20000])
    %     title('T21','FontName','Times','FontSize',14)
    %     subplot(1,3,3); imagesc(time,F,T22); 
    %     colormap(black); set(gca,'YDir','normal'); ylim([0 20000])
    %     title('T22','FontName','Times','FontSize',14)
    %     %Plotting - - - -

    % Computing angles and anisotropy measure:-----------------------------
    alphas = zeros(size(X)); C = zeros(size(X));
    for m = 1:size(X,2) %time
        for k = 1:size(X,1) %frequency
            T = [T11(k,m) T12(k,m); T21(k,m) T22(k,m)];
            [eig_vec, eig_val] = eig(T);

            % Angle of orientation:
            alphas(k,m) = atan(eig_vec(2,1)/eig_vec(1,1));

            % Anisotropy measure:
            if (eig_val(1) + eig_val(4)) > epsilon
                C(k,m) = ((eig_val(4) - eig_val(1)) / (eig_val(1) + eig_val(4)))^2;
            end            
        end        
    end

    R = fs^2/(hop*(size(X,1)-1)*2)*tan(alphas); 

    %   % Harmonic-Percussive-Residual Separation:---------------------------------
    %     rh = 20000; rp = 20000; c = 0.02;    
    %     % Computing masks:
    %     Mh = zeros(size(S)); Mp = Mh; 
    %     Mh(abs(R)<=rh&C>c) = 1; Mp(abs(R)>rp&C>c) = 1; Mr = ones(size(S))-Mh-Mp;
    %     % Computing STFTs of harmonic, percussive, and residual components:
    %     Xh = Mh.*X; Xp = Mp.*X; Xr = Mr.*X;

    %--------------------------------------------------------------------------

    
    alphasFChT = R./repmat(F(:),1,size(X,2));
    alphasFChT(1,:) = -10^10;
    alphasEST  = zeros(length(time), NumOfSources); 

    % Gaussian Kernel:
    kernel = @(x)(exp(-.5*x.^2)/sqrt(2*pi));

    for m = 1:size(X,2)
        % Defining PDF:
        if sum(C(:,m))~=0
            pdfest = @(x)((C(:,m))'*kernel((x-alphasFChT(:,m))/h)/h)/sum(C(:,m)); % Weighting by C 
        else
            pdfest = @(x)mean(kernel((x-alphasFChT(:,m))/h)/h);
        end

        % Points: (how many points?)
        %[xpdf,ypdf] = fplot(pdfest, [-5,5], 'k', 100);
        [xpdf,ypdf] = fplot(pdfest, [-10,10], 'k', 100);

%         if m == 15
%             plot_alphas_dist;
%             return
%         end
        
        % Estimating alpha
        if m == 1
            alphasEST(m,:) = zeros(1,NumOfSources);
        else
            alphasEST(m,:) = find_alphas(xpdf,ypdf, percent_trsh, NumOfSources);
        end        
    end

    alphasEST_orig = alphasEST;
    if filter_alphas_flag
        alphasEST = filter_alphas_median(alphasEST, ord_median);
        alphasEST = filter_alphas(alphasEST, ord);
    end

    alphasEST(abs(alphasEST)<10^-8) = 0;    

    % Comparing to theoretical values
    % Plotting-------------------------------------------------------------
% 
%     figure; plot(alphasEST_orig, 'xb', 'linewidth', 2);
%     axis([60 120 -1 1])
%     figure; plot(alphasEST, '+r', 'linewidth', 2);
%     axis([60 120 -1 1])
%     return
%     
    %     figure;
    %     plot(Temp,alphasTEO,'linewidth',1.5); hold on
    %     plot(Temp,alphasTEO2,'r','linewidth',1.5);
    %     plot(Temp,alphasEST(:,1),'ok','linewidth',1.5)
    %     plot(Temp,alphasEST(:,2),'ok','linewidth',1.5)
    %     hold off
    %     legend('Theoretical 1', 'Theoretical 2', 'Estimated 1', 'Estimated 2');
    % % Relative error:
    % figure; plot(Temp,(alphasTEO-alphasEST)./alphasTEO); hold on;
    % plot(Temp,alphasTEO,'k','linewidth',1.5);
    % plot(Temp,alphasEST,'r','linewidth',1.5)
    % hold off;
    % ylim([-5 5]); grid
    % legend('Ratio','Theoretical','Estimated')   

    % Computing theoretical values of alpha:    
    if useTheoAlphas
        if NumOfSources == 1            
            alphasTEO = computeAlphasTEO(x1,fs,f1(1,:), (size(X,1)-1)*2, hop);
        elseif NumOfSources == 2
            alphasTEO = [computeAlphasTEO(x1, fs, f1(1,:), Nf, hop)  computeAlphasTEO(x2, fs, f2(1,:), Nf, hop)];
        end
        %alphasEST = alphasTEO;
    end

    % Compute STFChTs using alphas: - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tic_STFT = tic;
    [X_2, F_2, time_2] = spectrogram(x, hanning(Nf_2), Nf_2 - hop_2, Nf_2, fs); 
    F_final = F_2;
    X_2 = abs(X_2);
    toc(tic_STFT)
    
    Num_frames_2 = length(time_2); 
    STFChT = zeros(Nf_2/2+1, Num_frames_2, NumOfSources);
    alphasEST_interp = zeros(length(time_2), NumOfSources);

    for alphas_ind = 1:NumOfSources
        % Interpolation of the estimated alphas   
        alphasEST_interp(:,alphas_ind) = interp1(time, alphasEST(:, alphas_ind), time_2, 'pchip');    
        for ii = 1:Num_frames_2
            t0 = (ii-1)*hop_2+1;
            xn = x(t0:t0+Nf_2-1);
            [~,mag_specs] = calculaFChT(xn, fs, Nf_2, 1, alphasEST_interp(ii, alphas_ind));
            STFChT(:,ii,alphas_ind) = abs(mag_specs(1:Nf_2/2+1));
        end
    end

    if combine_original_TFR
        [m_grid, k_grid] = meshgrid(time, F);
        [m_grid_2, k_grid_2] = meshgrid(time_2, F_2);
        X_1 = interp2(m_grid, k_grid, X_original, m_grid_2, k_grid_2, 'linear', 0);
        cat_spectrograms = cat(3, X_1, STFChT);
    else
        cat_spectrograms = STFChT;
    end

    if combine_large_STFT
        cat_spectrograms = cat(3, X_2, cat_spectrograms);
    end

    if combine_extra_large_STFT
        [X_3, F_3, time_3] = spectrogram(x, hanning(Nf_3), Nf_3 - hop_2, Nf_3, fs); 
        X_3 = abs(X_3);
        
        [m_grid, k_grid] = meshgrid(time_3, F_3);
        [m_grid_2, k_grid_2] = meshgrid(time_2, F_2);
        cat_spectrograms_interp = X_3;
        
        for ind = 1:size(cat_spectrograms,3)
            cat_spectrograms_interp = cat(3, cat_spectrograms_interp, interp2(m_grid_2, k_grid_2, ...
                                                        cat_spectrograms(:,:,ind), m_grid, k_grid, 'linear', 0));
        end
        F_final = F_3;
        cat_spectrograms = cat_spectrograms_interp;
    end
    
    cat_spectrograms = abs(cat_spectrograms( F_final <= ylim_vec(2)*1000, :, :)).^2;
%     cat_spectrograms = abs(cat_spectrograms( F_final <= ylim_vec(2)*1000, :, :));
    F_final = F_final(F_final <= ylim_vec(2)*1000);
    
    if NumOfSources > 1
        switch comb_method
            case 1 % SWGM
                STFChT_comb = SWGM_comb(cat_spectrograms, beta);
                % Saving TFR in a .mat file
                save([TFR_files_dir '/FChT_ST_LS/' signal_name], 'STFChT_comb', 'time_2', 'F_final');
                
            case 2 % LS
                STFChT_comb = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
                % Saving TFR in a .mat file
%                 save([TFR_files_dir '/FChT_ST_LS/' signal_name], 'STFChT_comb', 'time_2', 'F_final');                
                
            case 3 % SLS
                STFChT_comb = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
                % Saving TFR in a .mat file
%                 save([TFR_files_dir '/FChT_ST_SLS/' signal_name], 'STFChT_comb', 'time_2', 'F_final');

            case 4 % All methods
                
                STFT_1024 = cat_spectrograms(:,:,2);
                STFT_2048 = cat_spectrograms(:,:,1);
                
                % Saving STFTs in a .mat file
                save([path_check([TFR_files_dir '/STFT-1024/']) signal_name], 'STFT_1024', 'time_2', 'F_final');
                save([path_check([TFR_files_dir '/STFT-2048/']) signal_name], 'STFT_2048', 'time_2', 'F_final');

                % SWGM - - - - - - - - -
                STFChT_comb = SWGM_comb(cat_spectrograms, beta);
                
                % Saving TFR in a .mat file
                save([path_check([TFR_files_dir '/FChT_ST_SWGM/']) signal_name], 'STFChT_comb', 'time_2', 'F_final');
                
                % LS - - - - - - - - - 
                zeta = 200; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
                eta = 1; % controls the energy weighting process
                STFChT_comb = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
                
                % Saving TFR in a .mat file
                save([path_check([TFR_files_dir '/FChT_ST_LS/']) signal_name], 'STFChT_comb', 'time_2', 'F_final');
                
                % SLS - - - - - - - - - 
                zeta = 100; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
                eta = 1; % controls the energy weighting process
                STFChT_comb = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
                
                % Saving TFR in a .mat file
                save([path_check([TFR_files_dir '/FChT_ST_SLS/']) signal_name], 'STFChT_comb', 'time_2', 'F_final');

        end
    else
        STFChT_comb = STFChT;
    end
    
    toc(tic_all)
% end

% Plot_TFRs

% return


%% *********************************************************

%%%%%%%%%%%%%%%
%       Plotting   Thesis
%%%%%%%%%%%%%%%

% plot_spec('synth_vibrato_STFChT_2048', compress_dB_norm(cat_spectrograms(:,:,3), 80), time_2, F_final, plot_par, figs_path, []);
% 
% plot_spec('synth_vibrato_STFT_2048', compress_dB_norm(cat_spectrograms(:,:,1), 80), time_2, F_final, plot_par, figs_path, []);

plot_spec([signal_name '_STFT_2048'], compress_dB_norm(cat_spectrograms(F_final <= ylim_vec(2)*1000, :, 1), 80), ...
                                                                                    time_2, F_final, plot_par, figs_path, []);

plot_spec([signal_name 'Comb_TFR_FChTs_SLS_2048'], compress_dB_norm(STFChT_comb, 80), time_2, ...
                                                                            F_final(F_final <= ylim_vec(2)*1000), plot_par, figs_path, []);

return

% Comparing STFT vs. STFChT - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
range = 80; % range in dB for plotting

figure; subplot(1,2,1); imagesc(time_2, F_final/1000, ...
                              compress_dB_norm(cat_spectrograms(F_final <= ylim_vec(2)*1000, :, 1), range));

set(gca,'YDir','normal'); ylim(ylim_vec); colormap(abs(1-gray));
title('STFT')
xlabel('Time [s]')
ylabel('Frequency [kHz]')
pos = get(gca, 'Position');
pos(1) = .14; 
pos(3) = .41;
set(gca, 'Position', pos)

subplot(1, 2, 2); imagesc(time_2, F_final(F_final <= ylim_vec(2)*1000)/1000, ...
                                                                                            compress_dB_norm(STFChT_comb, range));
% subplot(1, 2, 2); imagesc(time_2, F_final(F_final <= ylim_vec(2)*1000)/1000, ...
%                                                                                             compress(STFChT_comb.^0.5));

set(gca,'YDir','normal'); ylim(ylim_vec); colormap(abs(1-gray));
if comb_method == 3
    title('SLS combination') % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    tit = ['TFRs_' signal_name '_Nsources_' num2str(NumOfSources) '_SLS'];
elseif comb_method == 2
    title('LS combination') % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    tit = ['TFRs_' signal_name '_Nsources_' num2str(NumOfSources) '_LS'];
end
xlabel('Time [s]')
set(gca,'YTick',[]);
pos = get(gca, 'Position');
pos(1) = .555;
pos(3) = .41;
set(gca, 'Position', pos)

% tit = ['TFRs_' signal_name '_Nsources_' num2str(NumOfSources)];

%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
% figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % D.Sc.
figProp = struct('size', 35, 'font','Times','lineWidth',2,'figDim',[1 1 900 500]); % D.Sc.
figFileName = [figs_path tit];
formatFig(gcf, figFileName, 'en', figProp);


return


%%%%%%%%%%%%%%%
%       Plotting   JAES
%%%%%%%%%%%%%%%

close all;

%     % Plots the FChT instances (don't think it is useful for the paper) - - - - - - - -
%     figure;
%     for ii = 1:NumOfSources
%         subplot(1,NumOfSources,ii); black = abs(1-gray); imagesc(time_2, F_2, compress(STFChT(:,:,ii)));
%         set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
%         title(sprintf('Short-Time Fan-Chirp Transform %d', ii),'FontName','Times','FontSize', 14)
%     end

% Comparing STFT vs. STFChT - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
figure; subplot(1,2,1); imagesc(time_2, F_final/1000, ...
                              compress(abs(cat_spectrograms(F_final <= ylim_vec(2)*1000, :, 1))));

set(gca,'YDir','normal'); ylim(ylim_vec); colormap(abs(1-gray));
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (kHz)')
pos = get(gca, 'Position');
pos(1) = .14; 
pos(3) = .41;
set(gca, 'Position', pos)

subplot(1, 2, 2); imagesc(time_2, F_final(F_final <= ylim_vec(2)*1000)/1000, ...
                                                                                            compress(STFChT_comb));

set(gca,'YDir','normal'); ylim(ylim_vec); colormap(abs(1-gray));
title('Combined TFR') % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
xlabel('Time (s)')
%     ylabel('Frequency (Hz)')
set(gca,'YTick',[]);
pos = get(gca, 'Position');
%     pos(1) = .58;
%     pos(3) = .41;
pos(1) = .555;
pos(3) = .41;
set(gca, 'Position', pos)

%     tit = ['TFRs_' signal_name '_Nsources_' num2str(NumOfSources) 'TEST'];
tit = ['TFRs_' signal_name '_Nsources_' num2str(NumOfSources)];

%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
% figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % D.Sc.
figProp = struct('size', 23,'font','Times','lineWidth',2,'figDim',[1 1 800 500]); % D.Sc.
figFileName = [figs_path tit];
formatFig(gcf, figFileName, 'en', figProp);

%     % Plots alphas  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%     figure;
%     plot(time, alphasEST_orig, 'xb');
%     if filter_alphas_flag
%         hold on;
%         plot(time_2, alphasEST_interp, '+r');
%         legend('Estimated', 'Filtered')
%     end    
%     if plotTheoAlphas
%         hold on;
%         plot(time, resMult.*alphasTEO, 'ok');
%         legend('Estimated', 'Theoretical')
%     end
%     hold off;
%     xlabel('Time (s)')
%     ylabel('\alpha')    
% %     tit = ['Alphas_' signal_name '_Nsources_' num2str(NumOfSources) 'TEST'];
%     tit = ['Alphas_' signal_name '_Nsources_' num2str(NumOfSources)];
%     
%     if isempty(ylim_alphas)
%         ylim_alphas = [-2,2];
%     end
%     
%     ylim(ylim_alphas)
%     
%     
%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 400]);
%     figFileName = [figs_path tit];
%     formatFig(gcf, figFileName, 'en', figProp);
