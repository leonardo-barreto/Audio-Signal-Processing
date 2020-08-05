% This function is intended to be a testing platform for the resampling algorithm. It performs sinusoidal analysis and calculates some pre-defined theoretical resampling curves,
% applying all of them to distort and then their inverted versions hoping to recover something as clos as possible to the original signal in the end.
%
%

function [outputResults] = MultipleCurveResampling(inputSignal,fs,windowType,windowSize,overlapPerc,fftPoints,filterCoeffs,waveCycles,alfa,DEBUG)

    [dummy,dummy2,inputSinAnalysisParameters] = SinusoidalAnalysis(inputSignal,fs,windowType,windowSize,overlapPerc,fftPoints,DEBUG);
    [inputDistortionCurves,dummy,outputResults.curveNames] = MakeDistortionCurves(inputSinAnalysisParameters,waveCycles,alfa);

    for curveIndex = 1:length(inputDistortionCurves)

        fprintf('\nThis is the ''%s'', curve no. %i of %i.\n',outputResults.curveNames{curveIndex},curveIndex, length(inputDistortionCurves));
        fprintf('Distorting...\n');

        outputResults.distortedSignals{curveIndex} = TimeVarying_Resample(inputSignal,inputSinAnalysisParameters,inputDistortionCurves{curveIndex},filterCoeffs,DEBUG);
        outputResults.distortionCurves{curveIndex} = inputDistortionCurves{curveIndex};

        [dummy,dummy2,distortedSinAnalysisParameters] = SinusoidalAnalysis(outputResults.distortedSignals{curveIndex},fs,windowType,windowSize,overlapPerc,fftPoints,DEBUG);
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