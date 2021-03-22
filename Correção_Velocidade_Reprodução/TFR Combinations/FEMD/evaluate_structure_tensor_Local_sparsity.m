clear variables; clc; close all;

%%%%%%%%%%%%%
%       Setting Paths       %
%%%%%%%%%%%%%

figs_path = path_check('./figs/'); 

if isunix
    addpath('./../../Audio/sinais_selecionados')
    addpath('./../../Audio')
    addpath('./Sinais')
%     addpath('./Resultados Saliencia')
else
    addpath('.\..\..\Audio')
    addpath('.\..\..\Audio\sinais_selecionados')
    addpath('.\Sinais')
%     addpath('./Resultados Saliencia')
end

%%%%%%%%%%%%%%%%%%
%       Parameters configuration     %
%%%%%%%%%%%%%%%%%%

% Analysis TFR - - - - - - - - - - - - - - - - - - - -
comb_method = 1; % 1- SWGM, 2 - LS, 3 - SLS

% - - - Plotting parameters - - -
ylim_vec = [0 20]; % in kHz
ylim_alphas = [];
black = abs(1-gray);

%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

%Get files' name (wav)
switch comb_method
    case 1 % SWGM
        TFR_files_dir = '/Volumes/HD_MAURICIO/datasets/MedleyDB/MedleyDB_TFRs/FChT_ST_SWGM_AES/';

    case 2 % LS
        TFR_files_dir = '/Volumes/HD_MAURICIO/datasets/MedleyDB/MedleyDB_TFRs/FChT_ST_LS_AES/';

    case 3 % SLS
        TFR_files_dir = '/Volumes/HD_MAURICIO/datasets/MedleyDB/MedleyDB_TFRs/FChT_ST_SLS_AES/';
end

TFR_file_name = files_in_path(TFR_files_dir,{'mat'});

for file_ind = 1:length(TFR_file_name)
    tic;
    %Load wav file - - - - - - - - - - - - - - - - - - - - -
    [data, fs_orig] = load(TFR_file_name(file_ind).name);
    
    signal_name = TFR_file_name(file_ind).clean_file_name;
    
    NumOfSources = 2;
    useTheoAlphas = 0; % For testing with controlled signals
    plotTheoAlphas = 0; % For testing with controlled signals
    
    % Make it mono
    if size(data,2) > 1
        x = mean(data.');
    else
        x = data;
    end
    x = x(:);
    
    % Resampling
    if fs_orig ~= fs
        x = resample(x,fs,fs_orig);
    end    
    

end

%% *********************************************************
%%%%%%%%%%
%       Plotting      %
%%%%%%%%%%
% close all;

%     % Plots the FChT instances (don't think it is useful for the paper) - - - - - - - -
%     figure;
%     for ii = 1:NumOfSources
%         subplot(1,NumOfSources,ii); black = abs(1-gray); imagesc(time_synth, F_synth, compress(STFChT(:,:,ii)));
%         set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
%         title(sprintf('Short-Time Fan-Chirp Transform %d', ii),'FontName','Times','FontSize', 14)
%     end


% Comparing STFT vs. STFChT - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
figure; subplot(1,2,1); imagesc(time_synth, F_synth(F_synth <= ylim_vec(2)*1000)/1000, ...
                                    compress(abs(S_synth(F_synth <= ylim_vec(2)*1000, :))));
%     figure; subplot(1,2,1); imagesc(time, F, compress(X_original)); % <<<<<<<<<<<<<<<<<<
set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
title('Spectrogram')
xlabel('Time (s)')
ylabel('Frequency (kHz)')
pos = get(gca, 'Position');
%     pos(1) = .14; 
%     pos(3) = .41;
pos(1) = .14; 
pos(3) = .41;
set(gca, 'Position', pos)

%     subplot(1,3, 2); imagesc(time_synth, F_synth, STFChT_comb);
%     set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
%     title('Short-Time Fan-Chirp Transform (SW)','FontName','Times','FontSize',14)
%     xlabel('Time (s)')
%     ylabel('Frequency (Hz)')

subplot(1, 2, 2); imagesc(time_synth, F_synth(F_synth <= ylim_vec(2)*1000)/1000, ...
                                                                                            compress(STFChT_comb));

set(gca,'YDir','normal'); ylim(ylim_vec); colormap(black);
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
figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 700 500]); % D.Sc.
figFileName = [figs_path tit];
formatFig(gcf, figFileName, 'en', figProp);

%     % Plots alphas  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%     figure;
%     plot(time, alphasEST_orig, 'xb');
%     if filter_alphas_flag
%         hold on;
%         plot(time_synth, alphasEST_interp, '+r');
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


%%     % Computing the entropy
S_synth = abs(S_synth(1:size(STFChT_comb,1),:));
S_interp = abs(S_interp(1:size(STFChT_comb,1),:));

%     sorted_samples = sort(S_interp(:));
sorted_samples = sort(S_synth(:));
sorted_samples = sorted_samples(end:-1:1);

cum_energy_vector = zeros(size(sorted_samples));
cum_energy_vector(1) = sorted_samples(1).^2;
for ind = 2:length(sorted_samples)    
    cum_energy_vector(ind) = cum_energy_vector(ind - 1) + sorted_samples(ind).^2;
end
cum_energy_vector =  cum_energy_vector/sum(sum(sorted_samples.^2));

energy_lims = [90 95 99];
energy_perc_spec1 = zeros(size(energy_lims));
energy_perc_spec2 = zeros(size(energy_lims));
energy_perc_comb = zeros(size(energy_lims));
GiniIndex = zeros(size(energy_lims,2), 3);

ind = 1;
for energy_lim = energy_lims/100   

    energy_vector_lim = cum_energy_vector(cum_energy_vector >= energy_lim);
    energy_sample = energy_vector_lim(1);
    lim = sorted_samples(cum_energy_vector == energy_sample);

%         bins = S_interp > lim;
    bins = S_synth > lim;

%         disp('Renyi:'); renyiEntropy = computeRenyiEntropy(cat(3, S_original(bins), STFChT_comb(bins)), ...
%                                                                                                             Nf, hop, renyi_ord)

    energy_perc_spec1(ind) = sum(sum(S_interp(bins).^2))/sum(sum(S_interp.^2));
    energy_perc_spec2(ind) = sum(sum(S_synth(bins).^2))/sum(sum(S_synth.^2));
    energy_perc_comb(ind) = sum(sum(STFChT_comb(bins).^2))/...
                                                        sum(sum(STFChT_comb.^2));

    GiniIndex(ind,:) = computeGiniIndex(cat(3, STFChT_comb(bins), S_interp(bins), S_synth(bins)));

    ind = ind + 1;        
%         disp('- - - - - - - - - - ')
end

GiniIndex(end,:)

figure; plot(energy_lims, 100*[(GiniIndex(:,1)./GiniIndex(:,2) - 1) (GiniIndex(:,1)./GiniIndex(:,3) - 1)]);
xlabel('Energy of selected bins for STFT-2048 (%)')
ylabel('Percentage gain in Gini Index')
legend('Comb TFR x STFT-1024', 'Comb TFR x STFT-2048')

tit = ['Gini_perc_gain_' signal_name];
%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 800 500]); % D.Sc.
figFileName = [figs_path tit];
formatFig(gcf, figFileName, 'en', figProp);

figure; plot(energy_lims, GiniIndex);
xlabel('Energy of selected bins for STFT-2048 (%)')
ylabel('Gini Index')
legend('Comb TFR', 'STFT-1024', 'STFT-2048')

tit = ['Gini_' signal_name];
%     figProp = struct('size', 25,'font','Times','lineWidth',2,'figDim',[1 1 800 900]); % IEEE
figProp = struct('size', 22,'font','Times','lineWidth',2,'figDim',[1 1 800 500]); % D.Sc.
figFileName = [figs_path tit];
formatFig(gcf, figFileName, 'en', figProp);

%     figure; plot(energy_lims, [energy_perc_spec1' energy_perc_spec2' energy_perc_comb']); 
%     xlabel('Energy of selected bins for STFT-2048 (%)')
%     ylabel('Energy of selected bins (%)')
%     legend('STFT-1024', 'STFT-2048', 'Combined TFR')
