function [frameArray,signalTrackArray,sinAnalysisParameters] = SinusoidalAnalysis(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints,DEBUG)

    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - STFT
    %
    %   3rd - Peak Detection and frequency Enhancement
    %
    %   4th - Construction of sinusoidal tracks
    %

    fprintf('\n\n------- SINUSOIDAL ANALYSIS STARTED ------\n\n');

    % -------------------------------------- STFT and spectrogram stage -------------------------------------------
        
        fprintf('Short-Time Fourier Transform starting...\n Sampling Rate(Hz): %i\n Window: %s (size %i, overlap %i%%) \n FFT Points: %i\n', samplingRate,windowType,windowSize,overlapPerc,fftPoints);

        [spectrgMatrix,freqComponents,timeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints);
        
        powerMatrixDB = 10*log10(powerMatrix); 

        totalFreqBins = length(freqComponents);
        totalFrames = length(timeInstants);

        % Building a signal frame as a peak detection entity
        signalFrame = {};
        signalFrame.totalFrames = totalFrames; %Total number of signal frames
        signalFrame.currentFrame = 1; %Current frame
        signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
        signalFrame.freqComponents = freqComponents; %frequency components vector
        

        %Outputting general parameters.
        sinAnalysisParameters.samplingRate = samplingRate;
        sinAnalysisParameters.timeInstants = timeInstants;
        sinAnalysisParameters.windowSize = windowSize;
        sinAnalysisParameters.totalFrames = totalFrames;
        sinAnalysisParameters.hopSize = floor(((100-overlapPerc)/100)*windowSize);

    % ---------------------------------------------- Peak Detection ------------------------------------------------

        fprintf('\nPeak Detection Starting...\n');

        detectedPeaksMatrix = {};
        detectedFrequenciesMatrix = {};

        if DEBUG == 1
            %Random frame chosen for DEBUG (temporary)
            DEBUG_FRAME = randi(totalFrames);
            %DEBUG_FRAME = 2;

            for frameCounter = 1:totalFrames
                signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
                signalFrame.currentFrame = frameCounter;
                frameArray(frameCounter) = signalFrame;
                if DEBUG_FRAME == signalFrame.currentFrame
                    [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,fftPoints,1);
                else
                    frameArray(frameCounter) = signalFrame;
                    [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,fftPoints,0);
                end
            end

        else

            for frameCounter = 1:totalFrames
                signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
                signalFrame.currentFrame = frameCounter;
                frameArray(frameCounter) = signalFrame;
                [detectedPeaksMatrix{frameCounter},detectedFrequenciesMatrix{frameCounter}] = DetectSpectralPeaks(signalFrame,samplingRate,fftPoints,0);
            end

        end

        for index = 1:totalFrames
            frameArray(index).peakMatrix = [detectedPeaksMatrix{index};detectedFrequenciesMatrix{index}];
        end

     % ---------------------------------------- Sinusoidal Tracking -----------------------------------------------

        fprintf('\nSinusoidal tracking starting...\n');

        % Array that holds all of the signal's tracks.
        % This array starts with an empty track, and is extended as the signal demands new tracks, maintaining
        % the information of tracks that ended. Later, another part of the algorithm will gather the starting
        % and ending frames of each track and organize them.
        signalTrackArray = [];
        signalTrackArray = setNewTrack();
        

        for frameCounter = 1:totalFrames

            signalTrackArray = PartialTracking2020(frameArray(frameCounter),signalTrackArray,2);

        end

        if DEBUG == 1
            organizedTracks = PlotTracks(frameArray,sinAnalysisParameters,signalTrackArray);
        end
    
    fprintf('\n\n------- SINUSOIDAL ANALYSIS FINISHED ------\n\n');

end