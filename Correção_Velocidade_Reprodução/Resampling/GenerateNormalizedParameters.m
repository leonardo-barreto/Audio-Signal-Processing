function [normalizedStepSize,numberNewSamples,normalizedCurrentPeriod] = GenerateNormalizedParameters(resampleFactors,windowIndex,originalPeriod,blockSize,DEBUG)

    %This repeats first and last factors.
    %This is used for the edges of the signal (to reach the first factor and to transition out of the last one.)
    %This is necessary because resampling blocks are defined in-between two subsequent window centers, and factors are associated with these centers.
    resampleFactorsAdjusted = [resampleFactors(1) resampleFactors resampleFactors(end)];

    currentPeriod = originalPeriod/resampleFactorsAdjusted(windowIndex);
    nextPeriod = originalPeriod/resampleFactorsAdjusted(windowIndex+1);

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
    b = 3*currentPeriod - nextPeriod - 2*blockSize*originalPeriod;
    c = -2*blockSize*originalPeriod;

    numberNewSamples = floor(max(roots([a b c])));
    stepSize = (nextPeriod - currentPeriod)/(numberNewSamples + 1);
    
    normalizedStepSize = stepSize/originalPeriod;
    normalizedCurrentPeriod = currentPeriod/originalPeriod;

end

