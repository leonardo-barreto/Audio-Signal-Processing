function enhancedFrequencies = EnhanceFrequency_DFT1(initialFrequencies,samplingRate)

    % This function computes a frequency enhancement for a detected peak by means of the DFT1 method.

    enhancedFrequencies = (samplingRate/pi)*sin(pi*initialFrequencies./samplingRate);

end