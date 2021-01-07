% This function is intended to be a testing platform for the resampling algorithm. It performs sinusoidal analysis and calculates some pre-defined theoretical resampling curves,
% applying all of them to distort and then their inverted versions hoping to recover something as clos as possible to the original signal in the end.
%
%

function [outputResults] = MultipleCurveResampling(inputSignal,samplingRate,windowSize,filterCoeffs,waveCycles,alfa,DEBUG)

    [dummy,dummy2,timeInstants,dummy3] = ComputeSTFT(inputSignal,samplingRate,'rectangular',windowSize,0,windowSize);

    inputSinAnalysisParameters.signalSize = length(inputSignal);
    inputSinAnalysisParameters.samplingRate = samplingRate;
    inputSinAnalysisParameters.timeInstants = timeInstants;
    inputSinAnalysisParameters.windowSize = windowSize;
    inputSinAnalysisParameters.totalFrames = length(timeInstants);
    inputSinAnalysisParameters.hopSize = windowSize;

    [inputDistortionCurves,outputResults.curveNames] = MakeDistortionCurves(inputSinAnalysisParameters,waveCycles,alfa);

    for curveIndex = 1:length(inputDistortionCurves)

        fprintf('\nThis is the ''%s'', curve no. %i of %i.\n',outputResults.curveNames{curveIndex},curveIndex, length(inputDistortionCurves));
        fprintf('Distorting...\n');

        outputResults.distortionCurves{curveIndex} = inputDistortionCurves{curveIndex};
        [outputResults.distortedSignals{curveIndex},outputResults.blockSizes{curveIndex}] = TimeVarying_ResampleDEBUG(inputSignal,inputSinAnalysisParameters,inputDistortionCurves{curveIndex},filterCoeffs,[],DEBUG);
      

        fprintf('\nInverting the distortion...\n');

        inverseCurves{curveIndex} = 1./inputDistortionCurves{curveIndex};
        outputResults.inverseCurves{curveIndex} = inverseCurves{curveIndex};
        [outputResults.revertedSignals{curveIndex},dummy] = TimeVarying_ResampleDEBUG(outputResults.distortedSignals{curveIndex},inputSinAnalysisParameters,inverseCurves{curveIndex},filterCoeffs,outputResults.blockSizes{curveIndex},DEBUG);
        

    end

    load handel
    x = y;
    load laughter
    lolSize = min(length(y),length(x));
    sound(y(1:lolSize)/8+x(1:lolSize)/10,Fs)
    clearvars x y Fs
    
end