function output = ErrorMeasurements(inputSignal,referenceSignal)

    signalSize = min([length(inputSignal) length(referenceSignal)]);

    diffError = zeros(1,signalSize);
    diffError_abs = zeros(1,signalSize);
    relativeDiffError = zeros(1,signalSize);
    relativeDiffError_abs = zeros(1,signalSize);
    SNR = zeros(1,signalSize);

    inputSignal = inputSignal(1:signalSize);
    referenceSignal = referenceSignal(1:signalSize);

    diffError = inputSignal - referenceSignal;
    diffError_abs = abs(inputSignal - referenceSignal);

    relativeDiffError(1:signalSize) = diffError(1:signalSize)/referenceSignal(1:signalSize); 
    relativeDiffError_abs(1:signalSize) = diffError_abs(1:signalSize)/abs(referenceSignal(1:signalSize));
    
    relativeSquareDiffError(1:signalSize) = power(diffError(1:signalSize),2)/power(referenceSignal(1:signalSize),2);

    SNR(1:signalSize) = 10*log(ones(1,signalSize)/relativeSquareDiffError(1:signalSize));

    output = [diffError;diffError_abs;relativeDiffError;relativeDiffError_abs;SNR];

end

    