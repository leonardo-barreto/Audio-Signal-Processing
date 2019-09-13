function PlotSpectrogram(magSpecData, param)
% Calculate STFT
% Input:
%   - wavData           : wave form mono
%   - para
%       - frameLen      : frame length in points
%       - hop           : hop size in points
%       - window        : a window function vector of frameLength
%       - fftLen        : fft length in points
% Output:
%   - SpecData          : complex spectrogram, dimension = [fftLength/2+1, total frame number]
%
% Author: Zhiyao Duan
% Last modified: 5/24/2013

fftLen = param.fftLen;
fs = param.fs;
nFrames = param.nFrames;

fAxis = (0:fftLen/2)*fs/fftLen;     % frequency axis
figure;
imagesc(1:nFrames, fAxis, 20*log10(magSpecData));
set(gca,'YDir','normal');
title('Spectrogram');
xlabel('Time (Frame number)');
ylabel('Frequency (Hz)');
% ylim([0 8000]);