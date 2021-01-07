function [distortionCurves,curveNames] = MakeDistortionCurves(sinAnalysisParameters,waveCycles,alfa)

    signalSize = sinAnalysisParameters.signalSize;
    samplingRate = sinAnalysisParameters.samplingRate;              
    timeInstants = sinAnalysisParameters.timeInstants;                      
    frameSize = sinAnalysisParameters.windowSize;                           
    totalFrames = sinAnalysisParameters.totalFrames;                       
    hopSize = sinAnalysisParameters.hopSize;

    
    signalDuration = signalSize/samplingRate;


    sineCurve.main = 1 + alfa*cos(2*pi*(waveCycles/signalDuration)*timeInstants - pi/2);

    triangCurve.main = 1 + alfa*sawtooth(2*pi*(waveCycles/signalDuration)*timeInstants + pi/2,1/2);

    %stepCurve.main = 1.5 - alfa*square(2*pi*(waveCycles/signalDuration)*timeInstants);
    stepCurve.main = (1+alfa) - alfa*rectpuls(timeInstants,max(timeInstants));


    distortionCurves = {sineCurve.main;triangCurve.main;stepCurve.main};

    curveNames = {'sineCurve' 'triangCurve' 'stepCurve'};

end