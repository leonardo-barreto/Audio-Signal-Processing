function [normalizedStepSize,numberOfSteps,normalizedCurrentPeriod,normalizedNextPeriod] = GenerateNormalizedParameters(resampleFactors,blockIndex,originalPeriod,blockSize)

    blockLength = (blockSize-1)*originalPeriod; %tempo do bloco de influência de cada fator.

    resampleFactorsAdjusted = [1,resampleFactors,resampleFactors(length(resampleFactors))];


    currentPeriod = originalPeriod/resampleFactorsAdjusted(blockIndex);

    if (blockIndex < length(resampleFactors))
        nextPeriod = originalPeriod/resampleFactorsAdjusted(blockIndex+1);
    else
        nextPeriod = currentPeriod;
    end

    a = currentPeriod + nextPeriod;
    b = 3*currentPeriod - nextPeriod - 2*blockSize*originalPeriod + 2*originalPeriod;
    c = 2*originalPeriod(1 - blockSize);
    
    numberOfSteps = floor((-b + sqrt(b^2 - 4*a*c))/2*a);

    stepSize = (nextPeriod - currentPeriod)/(numberOfSteps + 1);
    
    normalizedStepSize = stepSize/originalPeriod;
    normalizedCurrentPeriod = 1/resampleFactorsAdjusted(blockIndex);
    normalizedNextPeriod = 1/resampleFactorsAdjusted(blockIndex + 1);

end

