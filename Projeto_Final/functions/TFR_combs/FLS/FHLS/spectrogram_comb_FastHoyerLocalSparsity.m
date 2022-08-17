function [final_TFR, spectrograms_tensor_weighted] = spectrogram_comb_FastHoyerLocalSparsity(spectrograms_tensor, size_W_m_k, eta)
% spectrogram_comb is a function intended to combine different
% time-frequency representations by a weighted geometric mean method.
% V2: Ignore nulled time frames and do not compute the LS method
% Inputs:
%   - spectrograms_tensor: a tensor containing the aligned spectrograms stacked
%   in the third dimension
%   - size_w_S(1) and size_w_S(2): number (odd) of rows and columns of the 2D window

epsilon = 1e-10; % To avoid null sparsity

Hoyer_local_sparsity_measure_tensor = zeros(size(spectrograms_tensor));
LocalSparsity_ratio_all = zeros(size(spectrograms_tensor));

% Generate the 2D analysis window - - - -
if mod(size_W_m_k(1), 2) == 0    
    size_W_m_k(1) = size_W_m_k(1) + 1;
    fprintf('WARNING: size_w_S(1) must be an odd number! Using %d instead!\n', size_W_m_k(1));
end

if mod(size_W_m_k(2),2) == 0
    size_W_m_k(2) = size_W_m_k(2) + 1;
    fprintf('WARNING: size_w_S(2) must be an odd number! Using %d instead!\n', size_W_m_k(2));
end

% Symmetric hamming window
wr = window(@hamming, size_W_m_k(1));
wc = window(@hamming, size_W_m_k(2));
[maskr, maskc] = meshgrid(wc,wr);
window_S = maskc.*maskr;

for spec_ind = 1:size(spectrograms_tensor, 3)
    
    % Compute Local Energy
    local_energy = xcorr2(spectrograms_tensor(:,:,spec_ind), window_S) + epsilon;
    
    % Trim to remove borders
    local_energy = local_energy((size(window_S,1) - 1)/2 + 1 : end - (size(window_S,1) - 1)/2, ...
                                (size(window_S,2) - 1)/2 + 1 : end - (size(window_S,2) - 1)/2);

    % Compute Local Energy_2
    local_energy_2 = xcorr2(spectrograms_tensor(:,:,spec_ind).^2, window_S.^2) + epsilon;

    % Trim to remove borders
    local_energy_2 = local_energy_2((size(window_S,1) - 1)/2 + 1 : end - (size(window_S,1) - 1)/2, ...
                                (size(window_S,2) - 1)/2 + 1 : end - (size(window_S,2) - 1)/2);

    % Compute Hoyer Local Sparsity
    N = size_W_m_k(1)*size_W_m_k(2);
    Hoyer_local_sparsity_measure_tensor(:, :, spec_ind) = (sqrt(N) - (local_energy./sqrt(local_energy_2))) ./ ...
                 ((sqrt(N) - 1) .* local_energy.^.5);
             
    clear local_energy local_energy_2

end

% Combining spectrograms with weighted mean - - - - - - - - - - - - - - - - - - - -
disp('Combining spectrograms...')

% Sample Weighted mean
for spec_ind = 1:size(spectrograms_tensor,3)
    LocalSparsity_spects_aux = Hoyer_local_sparsity_measure_tensor + epsilon;
    LocalSparsity_spects_aux(:,:,spec_ind) = [];
    LocalSparsity_ratio_all(:,:,spec_ind) = (Hoyer_local_sparsity_measure_tensor(:,:,spec_ind)./prod(LocalSparsity_spects_aux, 3)).^eta;
    LocalSparsity_ratio_all(isnan(LocalSparsity_ratio_all(:,:,spec_ind))) = max(max(LocalSparsity_ratio_all(:,:,spec_ind)));
end

clear LocalSparsity_spects_aux
clear LocalSparsity_ratio

spectrograms_tensor_weighted = spectrograms_tensor.*LocalSparsity_ratio_all./sum(LocalSparsity_ratio_all,3);

spectrograms_tensor_weighted(isnan(spectrograms_tensor_weighted)) = 0;

final_TFR = sum(spectrograms_tensor_weighted,3);

orig_energy = sum(sum(spectrograms_tensor(:,:,1)));
comb_energy = sum(sum(final_TFR));

final_TFR = final_TFR*orig_energy/comb_energy;

end