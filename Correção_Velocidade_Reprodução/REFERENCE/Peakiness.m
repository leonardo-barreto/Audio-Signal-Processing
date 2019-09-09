function pk = Peakiness (X, Nw, Ndft)
% Computation of the peakiness feature of a single frame or block according
% with "THE TONALNESS SPECTRUM: FEATURE-BASED ESTIMATION OF TONAL
% COMPONENTS, DAFx2013, Kraft and Zolzer and Lerch"
% inputs:
%       X:      magnitude spectrum of signal frame
%       Nw:     window function length
%       Ndft:   fft buffer size (for zero-paddind)
% output:
%       pk:     peakiness feature

M = round(2*Ndft/Nw);   % spectral main lobe width converted into bins

pk = (X(1:end-2*M,:)+X(2*M+1:end,:))./X(1+M:end-M,:); % peakiness feature calculation

pk = [ones(M,size(X,2))*10 ; pk ; ones(M,size(X,2))*10];  % compensating the borders' samples

end