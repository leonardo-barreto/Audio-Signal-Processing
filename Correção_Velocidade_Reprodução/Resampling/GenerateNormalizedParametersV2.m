function [normalizedStepSize,numberNewSamples,normalizedCurrentPeriod,normalizedNextPeriod] = GenerateNormalizedParameters(resampleFactors,windowIndex,originalPeriod,blockSize)

    resampleFactorsAdjusted = [1,resampleFactors,resampleFactors(length(resampleFactors))];
    %Adiciona 1 para interpolar com o 1o fator, e repete o último para não interpolar nada na segunda metade da última janela.
    currentPeriod = originalPeriod/resampleFactorsAdjusted(windowIndex);
    nextPeriod = originalPeriod/resampleFactorsAdjusted(windowIndex+1);

    a = currentPeriod + nextPeriod;
    b = 3*currentPeriod - nextPeriod - 2*blockSize + 2;
    c = 2*(1 - blockSize);

    delta = power(b,2) - 4*a*c;
    
    numberNewSamples = floor(-b + sqrt(delta)/2*a);

    stepSize = (nextPeriod - currentPeriod)/(numberNewSamples + 1);
    
    normalizedStepSize = stepSize/originalPeriod;
    normalizedCurrentPeriod = 1/resampleFactorsAdjusted(windowIndex);
    normalizedNextPeriod = 1/resampleFactorsAdjusted(windowIndex + 1);

end

