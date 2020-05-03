function [timePositions,signalLength] = GenSamplingCurve (totalBlocks,blockSize)

%This function will try to generate a custom sampling curve in order to generate arbitrarily sampled signals.

baseSamplingRate = 48000; %block calculations have to interpret the signal as uniformly sampled. 
                          %Therefore, window length in seconds is given by this rate.

%windowLength = windowSize/baseSamplingRate;
%windowCenter = windowSize/2;
%blockSize = windowSize;

signalLength = 0; %time length
timePositions = [];
blocks = [];

%samples = 0;
%while (samples < blockSize/2)
%    currentPeriod = 1/baseSamplingRate;
%    signalLength = signalLength + currentPeriod;
%    timePositions(end+1) = signalLength;
%    samples = samples + 1;
%end

% GENERATE THE ARRAY OF SAMPLING RATES!

samplingRateArray = baseSamplingRate*ones(totalBlocks,1);
randomArray = -baseSamplingRate/12 + (2*baseSamplingRate/12)*rand(totalBlocks,1);

samplingRateArray = samplingRateArray + randomArray;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for currentRate = samplingRateArray
    samples = 0;
    while (samples < blockSize)
        currentPeriod = 1/currentRate;
        signalLength = signalLength + currentPeriod;
        timePositions(end+1) = signalLength;
        samples = samples + 1;
    end
end

