%This doesn't work
if (PHASE_CORRECT == 1)
    phaseCorrection = blockSize+1 - (blockNewSamples+1)*currentPeriod - (blockNewSamples+1)*stepSize*(blockNewSamples/2)
    phaseCorrectionMY = (blockSize+1) - blockNewSamples*(currentPeriod + stepSize*(blockNewSamples-1)/2)/originalPeriod
    newSamplePosition = newSamplePosition - phaseCorrection
end

%Phase Correction
if PHASE_CORRECT == 1
    phaseCorrection = 0; %Factor that will correct the synchrony between processed blocks.
end