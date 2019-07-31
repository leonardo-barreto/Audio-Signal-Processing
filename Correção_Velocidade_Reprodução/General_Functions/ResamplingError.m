function ResamplingError (errorMessage,systemStatus)

    [windowIndex,numberOfBlocks,stepSize,currentStep,numberOfSteps,newSamplePosition,inputSignalSize,outputSize,convolutionIndex_left,convolutionIndex_right] = systemStatus;

    disp(sprintf('Window %i of %i', windowIndex, numberOfBlocks)) 
    disp(sprintf('Step Size = %d', stepSize))
    disp(sprintf('Current Step: %i of %i',currentStep,numberOfSteps))
    format shortG
    disp(sprintf('Sample Position: %d', newSamplePosition))
    disp(sprintf('Original Size: %i', inputSignalSize))
    disp(sprintf('Resampled Size (up to now): %i', outputSize))
    disp(sprintf('Left Index = %i', convolutionIndex_left))
    disp(sprintf('Right Index = %i', convolutionIndex_right))
    error(errorMessage);

end