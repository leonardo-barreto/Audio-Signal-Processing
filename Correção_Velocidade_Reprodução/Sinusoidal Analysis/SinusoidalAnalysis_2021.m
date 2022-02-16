function [frameArray,signalTrackArray,TFParams] = TFAnalysis(signal_name)
    
    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - TFR
    %
    %   3rd - Peak Detection and frequency Enhancement
    %
    %   4th - Construction of sinusoidal tracks
    %

    DEBUG = 1;

    fprintf('\n\n------- SINUSOIDAL ANALYSIS STARTED ------\n\n');

    % -------------------------------------- TFR stage -------------------------------------------

        % MRFCI
            main_MRFCI_test;
                
            %Gathering TFRs from Mauricio
            inputSignal = x;
            samplingRate = fs;
            powerMatrix = TFR_comb;
            timeInstants = Time_combined;
            freqComponents = Freq_combined;
            %fftPoints = 2*length(Freq_combined);
            fftPoints = max(Nf);
            hopSize = hop;
            
        %TFR Information    
        powerMatrix = compress_dB_norm(powerMatrix, plot_range);
        powerMatrixDB = 10*log10(powerMatrix);
        %powerMatrixDB = powerMatrix;
        totalFreqBins = length(freqComponents);
        totalFrames = length(timeInstants);
    

        % Building a signal frame as a peak detection entity
        signalFrame = {};
        signalFrame.totalFrames = totalFrames; %Total number of signal frames
        signalFrame.currentFrame = 1; %Current frame
        signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
        signalFrame.freqComponents = freqComponents'; %frequency components vector
        

        %Outputting general parameters.
        TFParams = {};
        TFParams.signalSize = length(inputSignal);
        TFParams.samplingRate = samplingRate;
        TFParams.timeInstants = timeInstants;
        TFParams.windowSize = fftPoints;
        TFParams.totalFrames = totalFrames;
        TFParams.hopSize = hopSize;
        %TFParams.hopSize = floor(((100-overlapPerc)/100)*windowSize);

        
        fprintf(' Total frames: %i\n',TFParams.totalFrames);
        fprintf(' Number of frequency bins: %i\n',signalFrame.totalFreqBins);
        fprintf('\nShort-Time Fourier Transform Finished.\n');

        %if DEBUG == 1 %call plot
            %PlotSpectrogram(freqComponents,timeInstants,powerMatrixDB);
        %end

    % ---------------------------------------------- Peak Detection ------------------------------------------------

        fprintf('\nPeak Detection Starting...\n');

        

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrumDB = powerMatrixDB(:,frameCounter);
            signalFrame.powerSpectrumThreshold = [];
            signalFrame.peakMatrix = [];
            signalFrame.currentFrame = frameCounter;
            [signalFrame.peakMatrix,signalFrame.powerSpectrumThreshold] = DetectSpectralPeaks(signalFrame,samplingRate,fftPoints,DEBUG);
            frameArray(frameCounter) = signalFrame;
        end

        if DEBUG == 1
            availableFrames = [];
            for frameCounter = 1:totalFrames
                if (~isempty(frameArray(frameCounter).peakMatrix))
                    availableFrames(end+1) = frameCounter;
                end
            end

            %Random frame chosen for DEBUG
            DEBUG_FRAME = availableFrames(randi(length(availableFrames)));
            %DEBUG_FRAME = 350;
            PlotPeakDetection(TFParams,frameArray(DEBUG_FRAME));
        end

        fprintf('Peak Detection Finished.\n');

    % ---------------------------------------- Sinusoidal Tracking -----------------------------------------------

        fprintf('\nSinusoidal tracking starting...\n');

        % Array that holds all of the signal's tracks.
        % This array starts with an empty track, and is extended as the signal demands new tracks, maintaining
        % the information of tracks that ended. Later, another part of the algorithm will gather the starting
        % and ending frames of each track and organize them.

        signalTrackArray = [];
        signalTrackArray = setNewTrack();
        

        for frameCounter = 1:totalFrames

            signalTrackArray = PartialTracking(frameArray(frameCounter),signalTrackArray,2);

        end

        if DEBUG == 1
            organizedTracks = PlotTracks(frameArray,TFParams,signalTrackArray);
        end
    
    fprintf('\n\n------- SINUSOIDAL ANALYSIS FINISHED ------\n\n');
