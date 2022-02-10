function spec_prod_w = spectrogram_comb(spects_matrix, Beta)
% spectrogram_comb is a function intended to combine different
% time-frequency representations by a weighted geometric mean method.
% Inputs:
%   - spects_matrix: a matrix containing the aligned spectrograms stacked
%   in the third dimension
%   - Beta: parameter that controls the ammount of weight applied. Good
%   values should be 0<Beta<5.

% figure;
% imagesc(spects_matrix(1:300,:,1))
% set(gca,'ydir','normal');
% 
% figure;
% imagesc(spects_matrix(1:300,:,2))
% set(gca,'ydir','normal');

% Normalizing global energy
for spect_ind = 1:size(spects_matrix,3)
    norm_factor = (size(spects_matrix, 1)*size(spects_matrix, 2)) / sum(sum(spects_matrix(:, :, spect_ind).^2));
    spects_matrix(:,:, spect_ind) = spects_matrix(:, :, spect_ind)*norm_factor;
end

alphas = zeros(size(spects_matrix));

% Computing the alphas  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for spec_ind = 1:size(spects_matrix,3)
        
    % Performing geometric mean of the other spects_matrix
    spects_matrix_aux = spects_matrix;
    spects_matrix_aux(:,:,spec_ind) = [];
    spects_geomean = geomean(spects_matrix_aux,3);

    alphas(:,:,spec_ind) = (spects_geomean./spects_matrix(:,:,spec_ind)).^Beta;
    alphas(isnan(alphas(:,:,spec_ind))) = max(max(alphas(:,:,spec_ind)));
end


% Combining spects_matrix with wheited geometric average - - - - - - - - - - - - - - - - - - - -
spects_matrix_alphas = spects_matrix.^alphas;
spec_prod_w = prod(spects_matrix_alphas,3).^(1./sum(alphas,3));
spec_prod_w(isnan(spec_prod_w)) = 0;

orig_energy = sum(sum(spects_matrix(:,:,1).^2));
comb_energy = sum(sum(spec_prod_w.^2));

spec_prod_w = spec_prod_w*orig_energy/comb_energy;

% figure;
% imagesc(spec_prod_w(1:300,:));
% set(gca,'ydir','normal');

% figure;
% plot(spects_matrix(:,10,1).^.25);
% hold on;
% plot(spects_matrix(:,10,2).^.25);
% plot(spec_prod_w(:,10).^.25);
% hold off;

% figure;
% plot(spects_matrix(:,10,1));
% hold on;
% plot(spects_matrix(:,10,2));
% plot(spec_prod_w(:,10));
% hold off;
% 
% figure;
% plot(alphas(:,10,1));
% hold on;
% plot(alphas(:,10,2));
% hold off;

end