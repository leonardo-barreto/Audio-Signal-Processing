function [TFR_base,signalTrackArray,TFParams] = SinusoidalAnalysis(inputSignal,fs,TFR_method,varargin)
    
    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - TFR
    %
    %   3rd - Peak Detection and frequency Enhancement
    %
    %   4th - Construction of sinusoidal tracks
    %

    DEBUG = 1;

    fprintf('\n------- SINUSOIDAL ANALYSIS STARTED ------\n');

    % -------------------------------------- TFR stage -------------------------------------------
        fprintf('\nTime-frequency analysis starting...\n');
        if strcmp(TFR_method,'STFT')
            [TFR_base, freqComponents, timeInstants] = TFAnalysis_STFT(inputSignal,fs);
        elseif strcmp(TFR_method,'CQT')
            [TFR_base, freqComponents, timeInstants] = TFAnalysis_CQT(inputSignal,fs);
        elseif strcmp(TFR_method,'SWGM')
            [TFR_base, freqComponents, timeInstants] = TFAnalysis_SWGM(inputSignal,fs);
        elseif strcmp(TFR_method,'FLS')
            [TFR_base, freqComponents, timeInstants] = TFAnalysis_FLS(inputSignal,fs);
        elseif strcmp(TFR_method,'FEMD')
            [TFR_base, freqComponents, timeInstants] = TFAnalysis_FEMD(inputSignal,fs); 
        elseif strcmp(TFR_method,'MRFCI')
            [TFR_base, freqComponents, timeInstants] = TFAnalysis_MRFCI(inputSignal,fs);
        else
            error('Invalid TFR Method. Options are STFT, CQT, SWGM, FLS or MRFCI.')
        end
        
        % -------------------------- HPSS stage -------------------------------
        if ~isempty(varargin)
            if strcmp(varargin(1),'HPSS')
                fprintf('\nHarmonic-percussive separation starting...\n');
                HPSS_Params = varargin{2}
                [TFR_base, ~, ~] = Iterative_HPR_Separation(TFR_base, HPSS_Params.nFilterSS, HPSS_Params.nFilterTr,...
                                                            HPSS_Params.nIter, HPSS_Params.method, HPSS_Params.kernel_option);
                fprintf('Harmonic-percussive separation finished.\n');
            else
                error('Invalid optional 1st argument. Either empty or ''HPSS''.');
            end
        end
        
        %TFR Information    
        powerMatrix = 10*log10(TFR_base);
        %powerMatrix = TFR_base;
        totalFreqBins = length(freqComponents);
        totalFrames = length(timeInstants);
    

        % Building a signal frame as a peak detection entity
        signalFrame = {};
        %signalFrame.totalFrames = totalFrames; %Total number of signal frames
        signalFrame.currentFrame = 1; %Current frame
        signalFrame.totalFreqBins = totalFreqBins; %Total number of FFT bins
        signalFrame.freqComponents = freqComponents'; %frequency components vector
        
        %Outputting general parameters.
        TFParams = {};
        TFParams.freqComponents = freqComponents;
        TFParams.timeInstants = timeInstants;
        %TFParams.windowSize = NFFT;
        %TFParams.totalFrames = totalFrames;
        %TFParams.hopSize = hop;
        %TFParams.hopSize = floor(((100-overlapPerc)/100)*windowSize);

        fprintf(' Total frames: %i\n',totalFrames);
        fprintf(' Number of frequency bins: %i\n',signalFrame.totalFreqBins);
        fprintf(' Method: %s\n', TFR_method)
        fprintf('\nTime-frequency analysis finished.\n');

        %if DEBUG == 1 %call plot
            %PlotSpectrogram(freqComponents,timeInstants,powerMatrix);
        %end

    % ---------------------------------------------- Peak Detection ------------------------------------------------

        fprintf('\nPeak detection starting...\n');

        for frameCounter = 1:totalFrames
            signalFrame.currentFrame = frameCounter;
            signalFrame.powerSpectrum = powerMatrix(:,frameCounter);
            signalFrame.powerSpectrumThreshold = [];
            signalFrame.peakMatrix = [];
            [signalFrame.peakMatrix,signalFrame.powerSpectrumThreshold] = PeakDetection(signalFrame,fs,2*(totalFreqBins-1),DEBUG);
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
            %DEBUG_FRAME = availableFrames(randi(length(availableFrames)));
            DEBUG_FRAME = 500;
            PlotPeakDetection(frameArray(DEBUG_FRAME),totalFrames);
            title(sprintf('Quadro %i de %i (%s)',DEBUG_FRAME,totalFrames,TFR_method))
        end

        fprintf('Peak detection finished.\n');

    % ---------------------------------------- Sinusoidal Tracking -----------------------------------------------

        fprintf('\nSinusoidal tracking starting...\n');

        % Array that holds all of the signal's tracks.
        % This array starts with an empty track, and is extended as the signal demands new tracks, maintaining
        % the information of tracks that ended. Later, another part of the algorithm will gather the starting
        % and ending frames of each track and organize them.

        signalTrackArray = setNewTrack();
        signalTrackArray = signalTrackArray([]);

        for frameCounter = 1:totalFrames

            signalTrackArray = PartialTracking_2023(frameArray(frameCounter),totalFrames,signalTrackArray,0);

        end

        %if DEBUG == 1
            %organizedTracks = PlotTracks(frameArray,timeInstants);
        %end
    
    fprintf('\n------- SINUSOIDAL ANALYSIS FINISHED ------\n');
