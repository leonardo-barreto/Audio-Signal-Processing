function [TFR_base,TFParams,signalTracks] = SinusoidalAnalysis(inputSignal,fs,TFR_method,backwardsTracking,varargin)
    
    %   This function makes a full sinusoidal analysis of a given signal, using auxiliary functions for modularity.
    %
    %   1st - TFR
    %
    %   3rd - Peak Detection and frequency Enhancement
    %
    %   4th - Construction of sinusoidal tracks
    %

    DEBUG = 0;

    % -------------------------------------- TFR stage -------------------------------------------
        fprintf('Time-frequency analysis in progress...\n');
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
                fprintf('Harmonic-percussive separation in progress...\n');
                HPSS_Params = varargin{2}
                [TFR_base, ~, ~] = Iterative_HPR_Separation(TFR_base, HPSS_Params.nFilterSS, HPSS_Params.nFilterTr,...
                                                            HPSS_Params.nIter, HPSS_Params.method, HPSS_Params.kernel_option);
                fprintf('Harmonic-percussive separation done.\n');
            else
                error('Invalid optional 1st argument. Either empty or ''HPSS''.');
            end
        end
        
        %TFR Information    
        powerMatrix = 10*log10(TFR_base);
        %powerMatrix = TFR_base;
        frameSize = length(freqComponents);
        totalFrames = length(timeInstants);
        
        %Outputting general TF parameters.
        TFParams.method = TFR_method;
        TFParams.freqComponents = freqComponents;
        TFParams.timeInstants = timeInstants;

        fprintf(' Total frames: %i\n',totalFrames);
        fprintf(' Number of frequency bins: %i\n',length(freqComponents));
        fprintf(' Method: %s\n', TFR_method)
        fprintf('Time-frequency analysis done.\n');

    % ---------------------------------------------- Peak Detection ------------------------------------------------

        fprintf('Peak detection in progress...\n');

        for frameIdx = 1:totalFrames
            signalFrame.powerSpectrum = powerMatrix(:,frameIdx);
            signalFrame.powerSpectrumThreshold = [];
            signalFrame.peakMatrix = [];
            signalFrame.currentFrame = frameIdx;
            [signalFrame.peakMatrix,signalFrame.powerSpectrumThreshold] = PeakDetection(signalFrame,TFParams,fs,2*(frameSize-1),DEBUG);
            frameArray(frameIdx) = signalFrame;
        end

        if DEBUG == 1
            availableFrames = [];
            for frameIdx = 1:totalFrames
                if (~isempty(frameArray(frameIdx).peakMatrix))
                    availableFrames(end+1) = frameIdx;
                end
            end

            %Random frame chosen for DEBUG
            %DEBUG_FRAME = availableFrames(randi(length(availableFrames)));
            DEBUG_FRAME = 500;
            PlotPeakDetection(frameArray(DEBUG_FRAME),totalFrames);
            title(sprintf('Quadro %i de %i (%s)',DEBUG_FRAME,totalFrames,TFR_method))
        end

        fprintf('Peak detection done.\n');

    % ---------------------------------------- Sinusoidal Tracking -----------------------------------------------

        fprintf('Sinusoidal tracking in progress...\n');

        % Array that holds all of the signal's tracks.
        % This array starts with an empty track, and is extended as the signal demands new tracks, maintaining
        % the information of tracks that ended. Later, another part of the algorithm will gather the starting
        % and ending frames of each track and organize them.

        signalTracks = setNewTrack();
        signalTracks = signalTracks([]);

        if backwardsTracking == 1
            for frameIdx = totalFrames:-1:1
                signalTracks = PartialTracking_2023MQ(frameArray(frameIdx),TFParams,signalTracks,1);
            end
            for trackIdx = 1:length(signalTracks)
                signalTracks(trackIdx).powerEvolution = flip(signalTracks(trackIdx).powerEvolution);
                signalTracks(trackIdx).frequencyEvolution = flip(signalTracks(trackIdx).frequencyEvolution);
                signalTracks(trackIdx).currentPower = signalTracks(trackIdx).powerEvolution(end);
                signalTracks(trackIdx).currentFrequency = signalTracks(trackIdx).frequencyEvolution(end);
                trackFinalFrame =  totalFrames - (signalTracks(trackIdx).startFrame-1);
                signalTracks(trackIdx).startFrame = totalFrames - (signalTracks(trackIdx).finalFrame-1);
                signalTracks(trackIdx).finalFrame = trackFinalFrame;
            end
        else
            for frameIdx = 1:totalFrames
                signalTracks = PartialTracking_2023MQ(frameArray(frameIdx),TFParams,signalTracks,0);
            end
        end    
        fprintf('Sinusoidal tracking done.\n');
end
