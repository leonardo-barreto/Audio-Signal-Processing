function comb_spectrogram = RM_comb(spects_matrix)
% RM_comb is a function intended to combine different
% time-frequency representations by reciprocal mean method.
% Inputs:
%   - spects_matrix: a matrix containing the aligned spectrograms stacked
%   in the third dimension

comb_spectrogram = mean(spects_matrix.^(-1), 3).^(-1);

orig_energy = sum(sum(spects_matrix(:,:,1)));
comb_energy = sum(sum(comb_spectrogram));

comb_spectrogram = comb_spectrogram*orig_energy/comb_energy;

end