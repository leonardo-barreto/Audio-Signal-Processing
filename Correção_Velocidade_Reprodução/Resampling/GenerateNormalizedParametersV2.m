function [normalizedStepSize,numberNewSamples,normalizedCurrentPeriod] = GenerateNormalizedParametersV2(resampleFactors,windowIndex,originalPeriod,blockSize,DEBUG)

    %This adds 1 as a default first factor and repeats the last factor in the end.
    %This is used for the edges of the signal (to reach the first factor and to transition out of the last one.)
    %This is necessary because resampling blocks are defined in-between two subsequent window centers, and factors associated with these centers.
    resampleFactorsAdjusted(1) = 1;
    resampleFactorsAdjusted(2:length(resampleFactors)+1) = resampleFactors;
    resampleFactorsAdjusted(length(resampleFactorsAdjusted)+1) = resampleFactors(length(resampleFactors));

    currentPeriod = 1/resampleFactorsAdjusted(windowIndex);
    nextPeriod = 1/resampleFactorsAdjusted(windowIndex+1);

    if DEBUG == 1
        if (currentPeriod ~= nextPeriod)
            disp ('Transition!')
            X = sprintf('Current Period: %d', currentPeriod);
             disp(X)
            X = sprintf('Next Period: %d', nextPeriod);
            disp(X)
        end
    end

    a = currentPeriod + nextPeriod;
    b = 3*currentPeriod - nextPeriod - 2*blockSize + 2;
    c = 2*(1 - blockSize);

    numberNewSamples = floor(max(roots([a b c])));
    stepSize = (nextPeriod - currentPeriod)/(numberNewSamples + 1);
    
    normalizedStepSize = stepSize;
    normalizedCurrentPeriod = currentPeriod;

end

