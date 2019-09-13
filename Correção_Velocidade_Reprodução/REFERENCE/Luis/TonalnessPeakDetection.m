function [freqBuffer, magBuffer] = TonalnessPeakDetection (magSpecData, freqBuffer, param)

frameLen = param.frameLen;
fftLen = param.fftLen;
fs = param.fs;
nFrames = param.nFrames;
alpha_at = param.alpha_at;
specThreshFactor = param.specThreshFactor;
tonalnessTh = param.tonalnessTh;

% Tonalness computation
param = {};
param.frameLen = frameLen;
param.fftLen = fftLen;
param.alpha_at = alpha_at;

tonalnessMatrix = CalculateTonalness (magSpecData, param);
magBuffer = freqBuffer;
% Iterative peak detection
% for frameNum = 133:nFrames
for frameNum = 1:nFrames
    CurMagSpecData=magSpecData(:,frameNum);
    CurTonalness=tonalnessMatrix(:,frameNum);
    [~, locs] = findpeaks(CurMagSpecData,'MINPEAKHEIGHT',specThreshFactor*max(CurMagSpecData)); % local maxima above the threshold
    curPeaks = intersect(locs,find(CurTonalness > tonalnessTh)); % tonalness criteria
%     curPeaks = locs;
    curPeaks = curPeaks(:);       % converting into column vector;
%     PlotTonalness(CurTonalness,curPeaks,fftLen,fs,CurMagSpecData)
    [curInterpPeaks, curInterpMags] = PeakInterpolation (CurMagSpecData, curPeaks);
    [curInterpMagsSorted, ind] = sort(curInterpMags,'descend');
    curInterpPeaksSorted = curInterpPeaks(ind);
    curFrequencies = (curInterpPeaksSorted(1:min(60,length(curPeaks)))-1) * fs / fftLen;
    freqBuffer(1:length(curFrequencies),frameNum) = curFrequencies;
    magBuffer(1:length(curFrequencies),frameNum) = curInterpMagsSorted(1:min(60,length(curPeaks)));
end
idx1 = sum(freqBuffer~=0, 2)~=0;
freqBuffer = freqBuffer(idx1, :);

idx2 = sum(magBuffer~=0, 2)~=0;
magBuffer = magBuffer(idx2, :);

end