function final_TFR = spectrogram_comb_local_sparsity_V2(spects_matrix, size_w_S, size_w_E, zeta, eta)
% spectrogram_comb is a function intended to combine different
% time-frequency representations by a weighted geometric mean method.
% V2: Ignore nulled time frames and do not compute the LS method
% Inputs:
%   - spects_matrix: a matrix containing the aligned spectrograms stacked
%   in the third dimension
%   - size_w_S(1) and size_w_S(2): number (odd) of rows and columns of the 2D window

epslon = 10e-8; % To avoid null sparsity
size_w_S = ceil(size_w_S);
size_w_E = ceil(size_w_E);

if mod(size_w_S(1), 2) == 0    
    size_w_S(1) = size_w_S(1) + 1;
    fprintf('WARNING: size_w_S(1) must be an odd number! Using %d instead!\n', size_w_S(1));
end

if mod(size_w_S(2),2) == 0
    size_w_S(2) = size_w_S(2) + 1;
    fprintf('WARNING: size_w_S(2) must be an odd number! Using %d instead!\n', size_w_S(2));
end

if mod(size_w_E(1), 2) == 0
    size_w_E(1) = size_w_E(1) + 1;
    fprintf('WARNING: size_w_E(1) must be an odd number! Using %d instead!\n', size_w_E(1));
end

if mod(size_w_E(2),2) == 0
    size_w_E(2) = size_w_E(2) + 1;
    fprintf('WARNING: size_w_E(2) must be an odd number! Using %d instead!\n', size_w_E(2));
end

LocalSparsity_ratio_all = zeros(size(spects_matrix));
LocalSparsity =  zeros(size(spects_matrix)); % matrix that stores local sparsity
LocalEnergy =  zeros(size(spects_matrix)); % matrix that stores local energy
% LocalEnergy_aux =  zeros(size(spects_matrix,1), size(spects_matrix,2)); % matrix that stores local energy for chosen samples

% Generating the 2D windows
wr = window(@hamming, size_w_S(1));
wc = window(@hamming, size_w_S(2));
[maskr, maskc] = meshgrid(wc,wr);
window_S = maskc.*maskr;

wr = window(@hamming, size_w_E(1));
wc = window(@hamming, size_w_E(2));

wc(ceil(end/2) + 1 : end) = 0;
[maskr, maskc] = meshgrid(wc,wr);
window_E = maskc.*maskr;

% Zero-padding the spectrograms to properly apply the 2D windows at the edges
max_size_w = max(size_w_E, size_w_S);

zp_specs = cat(1, zeros((max_size_w(1)-1)/2,size(spects_matrix, 2), size(spects_matrix,3)), spects_matrix, ...
                   zeros((max_size_w(1)-1)/2,size(spects_matrix, 2), size(spects_matrix,3)));
zp_specs = cat(2, zeros(size(zp_specs, 1), (max_size_w(2)-1)/2, size(zp_specs,3)), zp_specs, ...
                   zeros(size(zp_specs, 1), (max_size_w(2)-1)/2, size(zp_specs,3)));
               
% Computing local sparsity - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Windowing spectrogram
Sparsity_vec_aux = zeros(1,size(spects_matrix, 3));

energy_frames = sum(spects_matrix(:,:,end), 1);

disp('Computing local sparsity and local energy...')

for col = (max_size_w(2)-1)/2 + 1 : size(spects_matrix, 2) + (max_size_w(2)-1)/2
    if energy_frames(col - (max_size_w(2)-1)/2) > 0  % presence of some signal
        for row = (max_size_w(1)-1)/2 + 1 : size(spects_matrix, 1) + (max_size_w(1)-1)/2
            for spec_ind = 1:size(spects_matrix, 3)
                regionMatrixS = zp_specs(row-(size_w_S(1)-1)/2:row+(size_w_S(1)-1)/2,...
                                                      col-(size_w_S(2)-1)/2:col+(size_w_S(2)-1)/2, spec_ind);
                windowedRegionSparsity = regionMatrixS.*window_S;

                regionMatrixE = zp_specs(row-(size_w_E(1)-1)/2:row+(size_w_E(1)-1)/2,...
                                                      col-(size_w_E(2)-1)/2:col+(size_w_E(2)-1)/2, spec_ind);
                windowedRegionEnergy = regionMatrixE.*window_E;

                % Computing local sparsity and energy
                Sparsity = computeGiniIndex(windowedRegionSparsity);
                LocalSparsity(row - (max_size_w(1)-1)/2, col - (max_size_w(2)-1)/2, spec_ind) = Sparsity;
                Sparsity_vec_aux(spec_ind) = Sparsity;

                LocalEnergy(row - (max_size_w(1)-1)/2, col - (max_size_w(2)-1)/2, spec_ind) = sum(sum(windowedRegionEnergy));
            end
        end
    end
end

LocalEnergy(LocalEnergy < epslon) = LocalEnergy(LocalEnergy < epslon) + epslon; % To avoid 0
LocalSparsity(LocalSparsity < epslon) = LocalSparsity(LocalSparsity < epslon) + epslon; % To avoid 0

LocalEnergy_spects_min = min(LocalEnergy, [], 3);

% Combining spectrograms with weighted mean - - - - - - - - - - - - - - - - - - - -
disp('Combining spectrograms...')

for spec_ind = 1:size(spects_matrix,3)
    LocalSparsity_spects_aux = LocalSparsity;
    LocalSparsity_spects_aux(:,:,spec_ind) = [];
    LocalSparsity_ratio_all(:,:,spec_ind) = (LocalSparsity(:,:,spec_ind)./prod(LocalSparsity_spects_aux, 3)).^zeta;
    LocalSparsity_ratio_all(isnan(LocalSparsity_ratio_all(:,:,spec_ind))) = max(max(LocalSparsity_ratio_all(:,:,spec_ind)));
end

clear LocalSparsity_spects_aux
clear LocalSparsity_ratio

for spec_ind = 1:size(spects_matrix,3)
    % Energy normalization
    LocalEnergy_ratio = LocalEnergy_spects_min./LocalEnergy(:,:,spec_ind);
    spects_matrix(:,:,spec_ind) = spects_matrix(:,:,spec_ind).*(LocalEnergy_ratio).^eta;
end

clear LocalEnergy_ratio

final_TFR = sum(spects_matrix.*LocalSparsity_ratio_all,3)./sum(LocalSparsity_ratio_all,3);
final_TFR(isnan(final_TFR)) = 0;

orig_energy = sum(sum(spects_matrix(:,:,1)));
comb_energy = sum(sum(final_TFR));

final_TFR = final_TFR*orig_energy/comb_energy;

end