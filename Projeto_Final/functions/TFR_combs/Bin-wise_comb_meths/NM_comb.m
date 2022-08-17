function comb_spectrogram = NM_comb(spects_matrix)
% NM_comb is a function intended to combine different
% time-frequency representations via numerical mean method.
% Input:
%   - spects_matrix: a matrix containing the aligned spectrograms stacked
%   in the third dimension

% Combining spectrograms using numerical mean
comb_spectrogram = mean(spects_matrix, 3);

end