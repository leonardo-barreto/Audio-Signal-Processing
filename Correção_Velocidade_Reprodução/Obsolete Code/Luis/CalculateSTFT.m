function [MagSpecData, nFrames] = CalculateSTFT(wavData, param, normFlag)
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

frameLen = param.frameLen;
hop = param.hop;
win = param.win;
fftLen = param.fftLen;

if nargin > 2
    if normFlag == 1
        win = win ./ frameLen;
    end
end

L = length(wavData);
nFrames = floor((L-frameLen)/hop + 1);            % the total number of frames

SpecData = zeros(fftLen/2+1, nFrames);            % to store the spectrum

for frameNum = 1:nFrames
%     fprintf('.');
    startp = hop*(frameNum-1) + 1;
    endp = startp + frameLen - 1;
    CurFrame = wavData(startp:endp) .* win;
    CurMagSpec = fft(CurFrame, fftLen);
    SpecData(:, frameNum) = CurMagSpec(1:fftLen/2+1);      % take only the first half spectrum
end

MagSpecData = abs (SpecData);
MagSpecData (MagSpecData<eps) = eps;           % switching zeros to epsilons to deal with logarithm