function [distortionCurves,inverseCurves,curveNames] = MakeDistortionCurves(sinAnalysisParameters,waveCycles,alfa)

    signalSize = sinAnalysisParameters.signalSize;
    samplingRate = sinAnalysisParameters.samplingRate;              
    timeInstants = sinAnalysisParameters.timeInstants;                      
    frameSize = sinAnalysisParameters.windowSize;                           
    totalFrames = sinAnalysisParameters.totalFrames;                       
    hopSize = sinAnalysisParameters.hopSize;

    
    signalDuration = signalSize/samplingRate;


    sineCurve.main = 1 + alfa*cos(2*pi*(waveCycles/signalDuration)*timeInstants);
    sineCurve.inverse = 2 - sineCurve.main;

    rampCurve.ascending = 1 + 2*alfa*timeInstants;
    rampCurve.descending = 1 - 2*alfa*timeInstants;

    triangCurve.main = 1 + 2*alfa*sawtooth(2*pi*(waveCycles/signalDuration)*timeInstants,1/2);
    triangCurve.main = triangCurve.main/triangCurve.main(1);
    triangCurve.inverse = 2 - triangCurve.main;

    distortionCurves = {sineCurve.main;rampCurve.ascending;triangCurve.main};
    inverseCurves = {sineCurve.inverse;rampCurve.descending;triangCurve.inverse};

    curveNames = {'sineCurve' 'rampCurve' 'triangCurve'};

end