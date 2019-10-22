function resampledSignal = CorrectVariationsMain(inputSignal,samplingRate)

    clc;
    fprintf('This is a program for corrections in the reprodution speed of audio signals.\n');
    fprintf('First, you have to enter the sinusoidal analysis parameters.\n');
    windowType = input('Window Type (hann or hamming): ');
    windowSize = input('Window Size (in samples, power of 2): ');
    overlapPerc = input('Overlap between windows (in %): ');
    fftPoints = input('Size of the FFT buffer (in samples, power of 2): ');
    DEBUG = input('Do you wish to receive debug messages (there will be a lot of them...)?(y/n) ');

    if strcmp(DEBUG,'y')
        DEBUG = 1;
    else
        DEBUG = 0;
    end

    frameArray = [];
    signalTrackArray = [];

    fprintf('\n\nNow, you need to enter the resampling parameters.\n');
    filterCoeffs = input('Number of sinc points for analog resampling (power of 2): ');
    fprintf('Would you like to offset the general pitch of the signal?\n');
    pitchOffset = input('If not, you should enter 0. If yes, enter by how much (%): ');

    fprintf('\nProcess will start. This should take a while (specially because of resampling).\n\n');

    outputSignal = [];

    [frameArray,signalTrackArray,sinAnalysisParameters] = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,DEBUG);
    
    if DEBUG == 1
        temp = PlotTracks(frameArray,sinAnalysisParameters,signalTrackArray);
    end
    
    resampleFactors = ExtractDeviationCurve (frameArray,signalTrackArray,DEBUG);
    timeRes = sinAnalysisParameters.timeInstants(10)-sinAnalysisParameters.timeInstants(9);
    trackTimes = 0:timeRes:sinAnalysisParameters.timeInstants(end)+2;
    figure;
    plot(trackTimes(1:sinAnalysisParameters.totalFrames),resampleFactors,'LineWidth',2);
    X = sprintf('Pitch Deviation Curve');
    title(X);
    xlabel('Time(s)');
    ylabel('Relative Frequency');

    resampleFactors = resampleFactors+(pitchOffset/100);
    
    resampledSignal = TimeVarying_Resample(inputSignal,sinAnalysisParameters,resampleFactors,filterCoeffs,DEBUG);

    fprintf('\n\nCorrections applied. Now, you should use the audiowrite function to create your audio file. Have fun!\n');

    load handel
    x = y;
    load laughter
    lolSize = min(length(y),length(x));
    sound(y(1:lolSize)/8+x(1:lolSize)/10,Fs)
    clearvars x y Fs

end

