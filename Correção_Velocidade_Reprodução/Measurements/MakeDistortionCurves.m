function [distortionCurves,inverseCurves,curveNames] = MakeDistortionCurves(sinAnalysisParameters,waveCycles,alfa)

    signalSize = sinAnalysisParameters.signalSize;
    samplingRate = sinAnalysisParameters.samplingRate;              
    timeInstants = sinAnalysisParameters.timeInstants;                      
    frameSize = sinAnalysisParameters.windowSize;                           
    totalFrames = sinAnalysisParameters.totalFrames;                       
    hopSize = sinAnalysisParameters.hopSize;

    
    signalDuration = signalSize/samplingRate;


    rampCurve.ascending = 1 + alfa*timeInstants;
    rampCurve.descending = 1 - 2*alfa*timeInstants;
    %rampCurve.descending = 2 - rampCurve.ascending;

    sineCurve.main = 1 + alfa*cos(2*pi*(waveCycles/signalDuration)*timeInstants);
    sineCurve.inverse = 2 - sineCurve.main;

    stepCurve.main = 1 + alfa*square(2*pi*(waveCycles/signalDuration)*timeInstants);
    stepCurve.inverse = 2 - stepCurve.main;

    triangCurve.main = 1 + 2*alfa*sawtooth(2*pi*(waveCycles/signalDuration)*timeInstants,1/2);
    triangCurve.main = triangCurve.main/triangCurve.main(1);
    triangCurve.inverse = 2 - triangCurve.main;

    distortionCurves = {rampCurve.ascending;sineCurve.main;stepCurve.main;triangCurve.main};
    inverseCurves = {rampCurve.descending;sineCurve.inverse;stepCurve.inverse;triangCurve.inverse};

    curveNames = {'rampCurve' 'sineCurve' 'stepCurve' 'triangCurve'};

end