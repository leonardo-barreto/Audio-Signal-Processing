function final_TFR = spectrogram_comb_FastHoyerLocalSparsity_cell(spectrograms_cell, final_time_scale, ...
                                                    final_frequency_scale, size_W_t_f, eta)
% spectrogram_comb is a function intended to combine different
% time-frequency representations by a weighted geometric mean method.
% V2: Ignore nulled time frames and do not compute the LS method
% Inputs:
%   - spectrograms_tensor: a tensor containing the aligned spectrograms stacked
%   in the third dimension
%   - size_w_S(1) and size_w_S(2): number (odd) of rows and columns of the 2D window

epsilon = 1e-10; % To avoid null sparsity

Hoyer_local_sparsity_measure_tensor = zeros(length(final_frequency_scale), length(final_time_scale), size(spectrograms_cell,1));
spectrograms_tensor = zeros(length(final_frequency_scale), length(final_time_scale), size(spectrograms_cell,1));
LocalSparsity_ratio_all = zeros(size(spectrograms_tensor));

for spec_ind = 1:size(spectrograms_cell, 1)
    
    % Generating the 2D windows - - - -
 
    hop = spectrograms_cell{spec_ind, 3}(2) - spectrograms_cell{spec_ind, 3}(1);
    freq_res = spectrograms_cell{spec_ind, 2}(2);

    size_W_m_k = ceil(size_W_t_f./[freq_res, hop]);
    
    if mod(size_W_m_k(1), 2) == 0    
        size_W_m_k(1) = size_W_m_k(1) + 1;
%         fprintf('WARNING: size_w_S(1) must be an odd number! Using %d instead!\n', size_W_m_k(1));
    end
    
    if mod(size_W_m_k(2),2) == 0
        size_W_m_k(2) = size_W_m_k(2) + 1;
%         fprintf('WARNING: size_w_S(2) must be an odd number! Using %d instead!\n', size_W_m_k(2));
    end    
    
    % Symmetric window
    wr = window(@hamming, size_W_m_k(1));
    wc = window(@hamming, size_W_m_k(2));
    [maskr, maskc] = meshgrid(wc,wr);
    window_S = maskc.*maskr;

%     % Asymmetric window
%     wr = window(@hamming, size_W_m_k(1));
%     wc = window(@hamming, size_W_m_k(2));
%     wc(ceil(end/2) + 1 : end) = 0;
%     [maskr, maskc] = meshgrid(wc,wr);
%     window_S = maskc.*maskr;
  
%     window_S = window_S / (sum(sum(window_S)));
%     window_S = window_S / (sum(sum(window_S))).^.5;
%     window_S = window_S * (size(window_S, 1) * size(window_S, 2)) / (sum(sum(window_S))).^.5;

%     sum(sum(window_S.^2))
    
    % Compute Local Energy
    local_energy = xcorr2(spectrograms_cell{spec_ind, 1}, window_S) + epsilon;
    
    % Trim to remove borders
    local_energy = local_energy((size(window_S,1) - 1)/2 + 1 : end - (size(window_S,1) - 1)/2, ...
                                (size(window_S,2) - 1)/2 + 1 : end - (size(window_S,2) - 1)/2);

% 	size(spectrograms_cell{spec_ind, 1})
%     size(local_energy)
%     figure; imagesc(local_energy)
%     pause; %close
    
    % Compute Local Energy_2
    local_energy_2 = xcorr2(spectrograms_cell{spec_ind, 1}.^2, window_S.^2) + epsilon;

    % Trim to remove borders
    local_energy_2 = local_energy_2((size(window_S,1) - 1)/2 + 1 : end - (size(window_S,1) - 1)/2, ...
                                (size(window_S,2) - 1)/2 + 1 : end - (size(window_S,2) - 1)/2);                        

%     figure; imagesc(local_energy_2)
%     pause; close

    % Interpolate local sparsity and spectrogram to the final tf scale
    [X,Y] = meshgrid(spectrograms_cell{spec_ind, 3}(:), spectrograms_cell{spec_ind, 2}(:));
    [Xq, Yq] = meshgrid(final_time_scale, final_frequency_scale);

    spectrograms_tensor(:, :, spec_ind) = interp2(X, Y, spectrograms_cell{spec_ind, 1}, Xq, Yq, 'linear', 0);
    local_energy = interp2(X, Y, local_energy, Xq, Yq, 'linear', 0);
    local_energy_2 = interp2(X, Y, local_energy_2, Xq, Yq, 'linear', 0);
    
    % Compute Hoyer Local Sparsity
    smallest_hop = spectrograms_cell{1, 3}(2) - spectrograms_cell{1, 3}(1);
    highest_freq_res = spectrograms_cell{end, 2}(2);

    size_W_m_k = ceil(size_W_t_f./[highest_freq_res, smallest_hop]);
    
    N = size_W_m_k(1)*size_W_m_k(2);
%     Hoyer_local_sparsity_measure_tensor(:, :, spec_ind) = (sqrt(N) - (local_energy./sqrt(local_energy_2))) ./ ...
%                  ((sqrt(N) - 1) .* (local_energy/N).^.5);
%     Hoyer_local_sparsity_measure_tensor(:, :, spec_ind) = (sqrt(N) - (local_energy./sqrt(local_energy_2))) ./ ...
%                  ((sqrt(N) - 1) .* local_energy.^.5);
    Hoyer_local_sparsity_measure_tensor(:, :, spec_ind) = (sqrt(N) - (local_energy./sqrt(local_energy_2))) ./ ((sqrt(N) - 1));
    
    clear X Y Xq Yq local_energy local_energy_2
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure; imagesc(Hoyer_local_sparsity_measure_tensor(:, :, spec_ind))
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     pause;   
end

% Combining spectrograms with weighted mean - - - - - - - - - - - - - - - - - - - -
disp('Combining spectrograms...')

% weight = Hoyer_local_sparsity_measure_tensor.^eta;
% 
% final_TFR = sum(spectrograms_tensor.*weight, 3)./sum(weight, 3);

for spec_ind = 1:size(spectrograms_tensor,3)
    LocalSparsity_spects_aux = Hoyer_local_sparsity_measure_tensor + epsilon;
    LocalSparsity_spects_aux(:,:,spec_ind) = [];
    LocalSparsity_ratio_all(:,:,spec_ind) = (Hoyer_local_sparsity_measure_tensor(:,:,spec_ind)./prod(LocalSparsity_spects_aux, 3)).^eta;
    LocalSparsity_ratio_all(isnan(LocalSparsity_ratio_all(:,:,spec_ind))) = max(max(LocalSparsity_ratio_all(:,:,spec_ind)));
end

clear LocalSparsity_spects_aux
clear LocalSparsity_ratio

final_TFR = sum(spectrograms_tensor.*LocalSparsity_ratio_all,3)./sum(LocalSparsity_ratio_all,3);
final_TFR(isnan(final_TFR)) = 0;

orig_energy = sum(sum(spectrograms_tensor(:,:,1)));
comb_energy = sum(sum(final_TFR));

final_TFR = final_TFR*orig_energy/comb_energy;
end