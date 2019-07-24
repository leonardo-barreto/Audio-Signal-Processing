function [diffError diffError_abs relativeDiffError relativeDiffError_abs SNR] = ErrorMeasurements(inputSignal,referenceSignal)

    
    signalSize = min([length(inputSignal) length(referenceSignal)]);
    inputSignal = inputSignal(1:signalSize);
    referenceSignal = referenceSignal(1:signalSize);

    diffError = inputSignal - referenceSignal;
    diffError_abs = abs(inputSignal - referenceSignal);

    relativeDiffError(1:signalSize) = diffError(1:signalSize)/referenceSignal(1:signalSize); 
    relativeDiffError_abs(1:signalSize) = diffError_abs(1:signalSize)/abs(referenceSignal(1:signalSize));
    
    relativeSquareDiffError(1:signalSize) = power(diffError(1:signalSize),2)/power(referenceSignal(1:signalSize),2);

    SNR(1:signalSize) = 10*log(ones(1,signalSize)/relativeSquareDiffError(1:signalSize));
end

    