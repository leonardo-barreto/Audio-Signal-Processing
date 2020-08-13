function tonalness = CalculateTonalness (X, param)
% Computation of the peakiness feature of a single frame or block according
% with "THE TONALNESS SPECTRUM: FEATURE-BASED ESTIMATION OF TONAL
% COMPONENTS, DAFx2013, Kraft and Zolzer and Lerch" 
% inputs:
%       X:      magnitude spectrum of signal frame
%       Nw:     window function length
%       Ndft:   fft buffer size (for zero-paddind)
% output:
%       pk:     peakiness feature

frameLen = param.frameLen;
fftLen = param.fftLen;
alpha_at = param.alpha_at;

vpk = Peakiness (X, frameLen, fftLen);

% eps_pk = log(2) / (mean(median(vpk))^2);

vat = AmpThreshold (X, alpha_at);

% eps_at = log(2) / (mean(median(vat))^2);

eps_at = 0.5298;
eps_pk = 0.1443;

tpk = exp(-1*eps_pk * vpk.^2);
tat = exp(-1*eps_at * vat.^2);

tonalness = tpk .* tat;
% tonalness = sts(vpk, epsilon_pk) .* sts(vat, epsilon_at);
% tonalness =  sts(vat, epsilon_at);



end