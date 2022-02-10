function comb_spectrogram = minimax_comb(spects_matrix)
% NM_comb is a function intended to combine different
% time-frequency representations via numerical mean method.
% Input:
%   - spects_matrix: a matrix containing the aligned spectrograms stacked
%   in the third dimension

% Combining spectrograms using numerical mean
comb_spectrogram = min(spects_matrix, [], 3);

orig_energy = sum(sum(spects_matrix(:,:,1)));
comb_energy = sum(sum(comb_spectrogram));

comb_spectrogram = comb_spectrogram*orig_energy/comb_energy;

end