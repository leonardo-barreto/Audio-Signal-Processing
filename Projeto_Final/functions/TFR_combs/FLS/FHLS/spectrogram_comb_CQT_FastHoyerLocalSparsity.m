function ok = spectrogram_comb_CQT_FastHoyerLocalSparsity(input_sig, fs, sig_name, ...
                                                                                bins_per_octave, fmax, fmin, hop, ...
                                                                                eta, size_w_S, ...
                                                                                plot_par, plot_flags)

% spectrogram_comb_CQT_LocalSparsity is a function that combines different resolution TFR
%   Detailed explanation goes here
% close all;

% size_w_S = ceil(size_w_S);
% 
% if mod(size_w_S(1), 2) == 0    
%     size_w_S(1) = size_w_S(1) + 1;
%     fprintf('WARNING: size_w_S(1) must be an odd number! Using %d instead!\n', size_w_S(1));
% end
% 
% if mod(size_w_S(2),2) == 0
%     size_w_S(2) = size_w_S(2) + 1;
%     fprintf('WARNING: size_w_S(2) must be an odd number! Using %d instead!\n', size_w_S(2));
% end
% 
% if isunix
%     dirbar = '/';
% else
%     dirbar = '\';
% end
% 
% figs_path = path_check(['.', dirbar, '..', dirbar, 'figs', dirbar,...
%                                     'Fast_Hoyer_Local_Sparsity', dirbar, sig_name, dirbar]);

                                
                                
                                
% addpath cqt_toolbox
% 
% % Computing CQT spectrograms - - - - - - - - - - - - - - - - - - - - - - 
% 
% % % Drop frequencies outside [fmin fmax]:
% % w1 = 1.1*(fmin/(fs/2));
% % w2 = 0.9*(fmax/(fs/2));
% % [B,A] = butter(6,[w1 w2]); 
% % input_sig = filtfilt(B,A,input_sig); 
% 
% % Computing rasterized complex coefficients - - - - - - - - - - - - - - -
% spectrograms_all = cell(length(bins_per_octave), 3);
% for ind = 1:length(bins_per_octave)
%     Xcqt = cqtPerfectRast(input_sig, fmin, fmax, bins_per_octave(ind), fs, 'q', 1, 'atomHopFactor',...
%                                       0.25, 'thresh', 0.0005, 'win', 'sqrt_blackmanharris');
% 
%     spec_CQT = abs(Xcqt.spCQT).^2;
% 
%     emptyHops = Xcqt.intParams.firstcenter/Xcqt.intParams.atomHOP;
%     maxDrop = emptyHops*2^(Xcqt.octaveNr - 1) - emptyHops;
%     droppedSamples = (maxDrop - 1)*Xcqt.intParams.atomHOP + Xcqt.intParams.firstcenter;
%     TimeVec = (1:size(spec_CQT, 2))*Xcqt.intParams.atomHOP - Xcqt.intParams.preZeros + droppedSamples;
% 
%     spec_CQT(:,(TimeVec < 0) | (TimeVec > length(input_sig))) = [];
%     TimeVec((TimeVec < 0) | (TimeVec > length(input_sig))) = [];
% 
%     freq_bins_CQT = Xcqt.fmin*2.^(((1:Xcqt.octaveNr*Xcqt.bins)-1)/Xcqt.bins);
% 
%     spectrograms_all{ind, 1} = spec_CQT;
%     spectrograms_all{ind, 2} = TimeVec;
%     spectrograms_all{ind, 3} = freq_bins_CQT;
% end
% 
% 
% % Interpolating CQTs - - - - - - - - - - - - - - - - - - - - - - - - - - 
% 
% finalTimeVec = 0:hop:TimeVec(end)/fs;
% 
% final_freq_bins_CQT = spectrograms_all{end, 3};
% 
% N_w_str = '';
% for CQT_ind = 1:length(bins_per_octave)
%     N_w_str = [N_w_str sprintf('_%d', bins_per_octave(CQT_ind))];
%     [X,Y] = meshgrid(spectrograms_all{CQT_ind, 2}/fs, spectrograms_all{CQT_ind, 3});
%     [Xq, Yq] = meshgrid(finalTimeVec, final_freq_bins_CQT);
%     
%     spectrograms_all{CQT_ind, 1} = interp2(X, Y, spectrograms_all{CQT_ind, 1}, Xq, Yq, 'linear', 0);
%     
%     if plot_flags(1)
%         plot_spec_CQT(sprintf('CQT_%dbins', bins_per_octave(CQT_ind)), compress_dB_norm(spectrograms_all{CQT_ind, 1},80), ...
%                                            finalTimeVec, final_freq_bins_CQT, Xcqt, plot_par, figs_path);
%     end
% end
% 
% 
% spectrograms = zeros(length(final_freq_bins_CQT), length(finalTimeVec), length(bins_per_octave));
% for ind = 1:length(bins_per_octave)
%     spectrograms(:,:,ind) = spectrograms_all{ind, 1};
% end
% 
% 
% % Normalizing CQTs - - - - - - - - - - - - - - - - - - - - - - - - - - -  
% 
% orig_energy = sum(input_sig.^2);
% 
% for ind = 1:size(spectrograms,3)
%     spec_energy = sum(sum(spectrograms(:,:,ind)));
%     spectrograms(:,:,ind) = spectrograms(:,:,ind)*orig_energy/spec_energy;
%     
% %     % Check energy
% %     sum(sum(spectrograms(:,:,ind)))
% end

[spectrograms, final_freq_bins_CQT, finalTimeVec, Xcqt] = CQT_spectrogram_tensor_prep(input_sig, fs, bins_per_octave, fmax, fmin, hop);

% Combining spectrograms with Fast Hoyer Local Sparsity method - - - - -

for ii_eta = 1:length(eta)
%         combined_TFR = spectrogram_comb_local_sparsity(spectrograms, size_w_S, size_w_E, zeta(ii_beta), eta(ii_gamma));
%         [combined_TFR, LocalEnergy_ratio, LocalSparsity] = spectrogram_comb_local_sparsity_2(spectrograms, size_w_S, size_w_E, zeta(ii_zeta), eta(ii_eta));
    combined_TFR_HLS = spectrogram_comb_FastHoyerLocalSparsity(spectrograms, size_w_S, eta);


%     fig_name = sprintf('FastHoyerLocalSparsity_CQT_eta%1.1f_NrS%d_NcS%d_%s', eta(ii_eta),...
%                                                                 size_w_S(1), size_w_S(2), N_w_str);
    fig_name = sprintf('FastHoyerLocalSparsity_CQT');

%     fig_name(fig_name == '.') = '_';
%         window_E_fig_name = 'LocalSparsity_DFT_window_E';
%         window_S_fig_name = 'LocalSparsity_DFT_window_S';
    if plot_flags(2)
        plot_spec_CQT(fig_name, compress_dB_norm(combined_TFR_HLS,80), finalTimeVec, final_freq_bins_CQT, Xcqt, plot_par, figs_path);
        title('FHLS')
% %             plot_window_S(size_w_S, window_S_fig_name, plot_par(end));
% %             plot_window_E(size_w_E, window_E_fig_name, plot_par(end));
    end
end

ok = 1;
end