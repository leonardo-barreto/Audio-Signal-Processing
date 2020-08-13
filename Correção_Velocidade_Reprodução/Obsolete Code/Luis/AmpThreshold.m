function at = AmpThreshold (X, alpha_at)
% Computation of the amplitude threshold feature of a single frame or block according
% with "THE TONALNESS SPECTRUM: FEATURE-BASED ESTIMATION OF TONAL
% COMPONENTS, DAFx2013, Kraft and Zolzer and Lerch" 
% inputs:
%       X:      magnitude spectrum of signal frame
%       alpha:  smoothing coefficient
% output:
%       at:     amplitude threshold feature

% b = 1 - alpha_at;      % numerator coefficients
% a = [1 -alpha_at];     % denominator coefficients

% applying the filter over the frequency in both the forward and backward
% direction to compensate for group delay
% at = filtfilt(b, a, X) ./ X;      % smoothed magnitude spectrum
% at = filtfilt(b, a, X);      % smoothed magnitude spectrum

thr = filtfilt(alpha_at, [1 -(1-alpha_at)], abs(X));
at = thr ./ (X + eps);
end
