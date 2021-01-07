% This function is intended to compose a testing platform for the resampling algorithm. It receives a reference signal and computes the difference between it and several resampled
% versions of it. The resampledSignals should be a cell vector of signals.
%

function [outputDifferences] = ComputeResamplingDifference(referenceSignal,resampledSignals,time,DEBUG)

    outputDifferences = {};

    

    for curveIndex = 1:length(resampledSignals)

        minLength = min([length(resampledSignals{curveIndex}),length(referenceSignal)]);

        outputDifferences{curveIndex}(1:minLength) = referenceSignal(1:minLength) - resampledSignals{curveIndex}(1:minLength);
        
        if DEBUG == 1
            figure;
            plot(time(1:minLength),abs(outputDifferences{curveIndex}));
            %X = sprintf('Diferenca curva %i',curveIndex);
            %title(X,'FontSize', 30);
            set(gca,'FontSize', 40);
        end

    end
    
    
end