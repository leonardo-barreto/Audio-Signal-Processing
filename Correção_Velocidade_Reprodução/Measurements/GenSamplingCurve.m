function [timePositions,samplingRateArray,signalLength] = GenSamplingCurve (baseSamplingRate,totalBlocks,blockSize)

%This function will try to generate a custom sampling curve in order to generate arbitrarily sampled signals.

%baseSamplingRate = 48000; %block calculations have to interpret the signal as uniformly sampled. 
                          %Therefore, window length in seconds is given by this rate.

signalLength = 0; %time length
timePositions = [];
blocks = [];

% GENERATE THE ARRAY OF SAMPLING RATES!

samplingRateArray = baseSamplingRate*ones(totalBlocks,1);

% Non-Uniform Function

%functionArray = 4000*sin(pi*(0:1/baseSamplingRate:(totalBlocks-1/baseSamplingRate)));
functionArray = 4000*sin(2*pi*(0:1/totalBlocks:1-1/totalBlocks));

samplingRateArray = samplingRateArray + functionArray;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for currentRate = transpose(samplingRateArray)
    samples = 0;
    while (samples < blockSize)
        currentPeriod = 1/currentRate;
        signalLength = signalLength + currentPeriod;
        timePositions(end+1) = signalLength;
        samples = samples + 1;
    end
end

