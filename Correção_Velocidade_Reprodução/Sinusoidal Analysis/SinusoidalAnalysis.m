function y = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,DEBUG)

    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - STFT
    %
    %   2nd - Frequency Enhancement
    %
    %   3rd - Peak Detection
    %
    %   4th - Tracking of sinusoidal components
    %


    % STFT and spectrogram stage

    fprintf('\nSinusoidal Analysis started.\n Sampling Rate(Hz): %i\n Window: %s (size %i, overlap %i%%) \n FFT Points: %i\n', samplingRate,windowType,windowSize,overlapPerc,fftPoints);

    [spectrgMatrix,freqComponents_cyclical,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);
    
    powerMatrixDB = 10*log(powerMatrix);
    freqComponents = (samplingRate/2*pi).*freqComponents_cyclical;

    totalFreqBins = length(freqComponents);
    totalFrames = length(timeInstants);

    % Building a signal frame as a peak detection entity
    signalFrame = {};
    signalFrame.totalFrames = totalFrames; %Total number of signal frames
    signalFrame.currentFrame = 1; %Current frame
    signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
    signalFrame.freqComponents = freqComponents; %frequency components vector

    % Peak Detection

    fprintf('\nPeak Detection Started.\n');

    detectedPeaksMatrix = zeros(totalFreqBins,totalFrames);

    if DEBUG == 1
        %Random frame chosen for DEBUG (temporary)
        DEBUG_FRAME = floor(rand(1,1)*(signalFrame.totalFrames-1) + 1);

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            if DEBUG_FRAME == signalFrame.currentFrame
                detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,1);
            else
                detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,0);
            end
        end

    else

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = (powerMatrixDB(:,frameCounter));
            signalFrame.currentFrame = frameCounter;
            detectedPeaksMatrix(:,frameCounter) = DetectSpectralPeaks(signalFrame,0);
        end

    end

    %TEMPORARY
    y = detectedPeaksMatrix;

    % Sinusoidal Tracking

    amplitudeMatrixDB = 10*log(abs(spectrgMatrix));
    MAXTRACKS = 100; %Maximum of tracks per frame.

    %Building what is a sinusoidal track
    sinusoidalTrack = {};
    sinusoidalTrack.amplitudeEvolution = zeros(totalFrames); %Contains the amplitude values of the track through its existence.
    sinusoidalTrack.frequencyEvolution = zeros(totalFrames); %Contains the frequency values of the track through its existence.
    sinusoidalTrack.startFrame = 0; %Starting frame for the track
    sinusoidalTrack.finalFrame = 0; %Ending frame for the track
    sinusoidalTrack.hysteresis = 0; %Hysteresis counter
    

    %Building a signal frame as a sinusoidal tracking entity
    signalFrame = {};
    signalFrame.totalFrames = totalFrames; %Total number of signal frames
    signalFrame.currentFrame = 1; %Current frame
    signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
    signalFrame.freqComponents = freqComponents; %frequency components vector

    signalFrame.sinusoidalTracks(1:MAXTRACKS) = sinusoidalTrack;


    sinusoidalTracksMatrix = zeros(totalFreqBins,totalFrames);
    
    for frameCounter = 1:totalFrames
        signalFrame.amplitudeSpectrumDB = (amplitudeMatrixDB(:,frameCounter));
        signalFrame.currentFrame = frameCounter;
        sinusoidalTracksMatrix(:,frameCounter) = DetectSinusoidalTracks(signalFrame,0);
    end

    
    

end