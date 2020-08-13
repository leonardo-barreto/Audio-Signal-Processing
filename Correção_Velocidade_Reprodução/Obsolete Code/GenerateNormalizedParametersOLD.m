function [normalizedStepSize,numberOfSteps,normalizedPeriod] = GenerateNormalizedParameters(resampleFactors,windowIndex,inputPeriod,blockSize,DEBUG)

    blockLength = blockSize*inputPeriod; %tempo do bloco de influência de cada fator.

    currentPeriod = inputPeriod/resampleFactors(windowIndex);

    if (windowIndex < length(resampleFactors))
        nextPeriod = inputPeriod/resampleFactors(windowIndex+1);
    else
        nextPeriod = inputPeriod/resampleFactors(windowIndex);
    end
    

    stepSize = ((nextPeriod^2) - (currentPeriod^2))/(2*blockLength-(currentPeriod+nextPeriod));

        if stepSize > 0
            numberOfSteps = floor(2*blockLength/(currentPeriod+nextPeriod));
        else
            numberOfSteps = ceil(2*blockLength/(currentPeriod+nextPeriod));
        end

    
    normalizedStepSize = stepSize/inputPeriod;
    normalizedPeriod = 1/resampleFactors(windowIndex);

end

