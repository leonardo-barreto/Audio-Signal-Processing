% This function is intended to be a testing platform for the resampling algorithm. It performs sinusoidal analysis and calculates some pre-defined theoretical resampling curves,
% applying all of them to distort and then their inverted versions hoping to recover something as clos as possible to the original signal in the end.
%
%

function [outputResults] = MultipleCurveResampling_BeforeChanges(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,filterCoeffs,waveCycles,alfa,DEBUG)

    [dummy,dummy2,timeInstants,dummy3] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);

    inputSinAnalysisParameters.signalSize = length(inputSignal);
    inputSinAnalysisParameters.samplingRate = samplingRate;
    inputSinAnalysisParameters.timeInstants = timeInstants;
    inputSinAnalysisParameters.windowSize = windowSize;
    inputSinAnalysisParameters.totalFrames = length(timeInstants);
    inputSinAnalysisParameters.hopSize = floor(((100-overlapPerc)/100)*windowSize);

    [inputDistortionCurves,dummy,outputResults.curveNames] = MakeDistortionCurves(inputSinAnalysisParameters,waveCycles,alfa);

    for curveIndex = 1:length(inputDistortionCurves)

        fprintf('\nThis is the ''%s'', curve no. %i of %i.\n',outputResults.curveNames{curveIndex},curveIndex, length(inputDistortionCurves));
        fprintf('Distorting...\n');

        outputResults.distortedSignals{curveIndex} = TimeVarying_Resample(inputSignal,inputSinAnalysisParameters,inputDistortionCurves{curveIndex},filterCoeffs,DEBUG);
        outputResults.distortionCurves{curveIndex} = inputDistortionCurves{curveIndex};

        

        distortedSinAnalysisParameters.signalSize = length(inputSignal);
        distortedSinAnalysisParameters.samplingRate = samplingRate;
        distortedSinAnalysisParameters.timeInstants = timeInstants;
        distortedSinAnalysisParameters.windowSize = windowSize;
        distortedSinAnalysisParameters.totalFrames = length(timeInstants);
        distortedSinAnalysisParameters.hopSize = floor(((100-overlapPerc)/100)*windowSize);

        [dummy,distortedInverseCurves,dummy2] = MakeDistortionCurves(distortedSinAnalysisParameters,waveCycles,alfa);

        fprintf('\nInverting the distortion...\n');

        outputResults.revertedSignals{curveIndex} = TimeVarying_Resample(outputResults.distortedSignals{curveIndex},distortedSinAnalysisParameters,distortedInverseCurves{curveIndex},filterCoeffs,DEBUG);
        outputResults.inverseCurves{curveIndex} = distortedInverseCurves{curveIndex};

    end

    load handel
    x = y;
    load laughter
    lolSize = min(length(y),length(x));
    sound(y(1:lolSize)/8+x(1:lolSize)/10,Fs)
    clearvars x y Fs
    
end