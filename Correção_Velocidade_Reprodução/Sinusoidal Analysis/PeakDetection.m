function y = PeakDetection(inputSignalWindow,samplingRate,windowSize,DEBUG)

    % This function aims to detect spectral peaks in a given window.
    %
    % Firstly, an initial all-peaks approach is taken.
    %
    % Secondly, a variable-length TPSW filtering builds a variable threshold.
    %
    % Finally, the frequency values are enhanced by use of the DFT1 method.
    %

    % In some moment, this will eventually happen
    enhancedPeaks = ComputeDFT1(peaksVector,samplingRate);
    
end