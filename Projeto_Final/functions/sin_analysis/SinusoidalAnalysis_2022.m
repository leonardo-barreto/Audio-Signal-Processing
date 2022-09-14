function [TFR_base,signalTrackArray,TFParams] = SinusoidalAnalysis(inputSignal,samplingRate,TFR_method)
    
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

        if strcmp(TFR_method,'STFT')

            %% - - - - - - Parameters configuration - - - - - - 
                hop = 128;
                N_w = 4096;
                NFFT = N_w;

            %% - - - - - - -  TFR computation - - - - - - -
                [spectrgMatrix, freqComponents, timeInstants] = spectrogram(inputSignal, hanning(N_w,'periodic'), N_w-hop, NFFT, samplingRate);
                TFR_base = power(abs(spectrgMatrix),2);%/NFFT;
            
        elseif strcmp(TFR_method,'FLS')

            %% - - - - - - Parameters configuration - - - - - -
                NFFT = 4096;
                hop_DFT = 0.0029; % Final time interval between frames in seconds (interp after all computations)
                N_w = [1, 2, 3, 4]*1024;
                hop = 128;
                overlap_short = N_w(1) - hop;
                % 2D Analysis window
                time_span = 0.03; % time span of the 2D window in seconds
                Ncols_DFT = ceil(time_span/hop_DFT);
                Nrows_DFT = 3*N_w(end)/N_w(1);
                size_W_m_k = [Nrows_DFT Ncols_DFT];
                % Contrast factor for weights
                gamma = 20;

            %% - - - - - - -  TFR computation - - - - - - -
                [spectrograms_tensor, freqComponents, timeInstants] = spectrogram_tensor_prep(inputSignal, samplingRate, N_w, NFFT, overlap_short);
                TFR_base = spectrogram_comb_FastHoyerLocalSparsity(spectrograms_tensor, size_W_m_k, gamma);

        elseif strcmp(TFR_method,'MRFCI')

            %% - - - - - - Parameters configuration - - - - - -
                % Analysis TFR - - - - - - - - - - - - - - - - - - - -
                hop = 128;
                block_size = 2; % Analysis block size in seconds

                % Synthesis TFR - - - - - - - - - - - - - - - - - - -
                N_alphas = 7; % Number of alphas (I, in the paper)
                Nf = [1024 2048 3072 4096]; % Window sizes
                NFFT = Nf(length(Nf));
                asym_flags = [0 0 0 1]; % Indicate which Nf will be used to compute Fan-chirp instances
                FChT_flags = [0 1 1 1]; % Indicate which Nf will be used to compute Fan-chirp instances

                % Structure tensor - - - - - - - - -
                C_limits = [0 .25 .75 1 ; 1 2 2 3]; % Defines the points of transition for the lambda weights
                range = 100; % range in dB for the ST

                % G parameters
                sigma_t = Nf(end)/(4*samplingRate); % in s
                sigma_f = 100;% in Hz

                Nf_structure_tensor = Nf(4); % To compute the structure tensor

            %% - - - - - - -  TFR computation - - - - - - -
                [TFR_base, STFT_tensor, timeInstants, freqComponents] = MRFCI(inputSignal, samplingRate, Nf, block_size,...
                FChT_flags, asym_flags, Nf_structure_tensor, hop, N_alphas, C_limits, range, sigma_t, sigma_f);
                freqComponents = freqComponents'; %Transposing because MRFCI function outputs 1 x N vector

        else
            error('Invalid TFR Method. Options are STFT, FLS or MRFCI.')
        end
            
        %TFR Information    
        powerMatrix = 10*log10(TFR_base);
        %powerMatrix = TFR_base;
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
        TFParams.freqComponents = freqComponents;
        TFParams.timeInstants = timeInstants;
        %TFParams.windowSize = NFFT;
        %TFParams.totalFrames = totalFrames;
        %TFParams.hopSize = hop;
        %TFParams.hopSize = floor(((100-overlapPerc)/100)*windowSize);

        
        fprintf(' Total frames: %i\n',totalFrames);
        fprintf(' Number of frequency bins: %i\n',signalFrame.totalFreqBins);
        fprintf('\nShort-Time Fourier Transform Finished.\n');

        %if DEBUG == 1 %call plot
            %PlotSpectrogram(freqComponents,timeInstants,powerMatrix);
        %end

    % ---------------------------------------------- Peak Detection ------------------------------------------------

        fprintf('\nPeak Detection Starting...\n');

        

        for frameCounter = 1:totalFrames
            signalFrame.powerSpectrum = powerMatrix(:,frameCounter);
            signalFrame.powerSpectrumThreshold = [];
            signalFrame.peakMatrix = [];
            signalFrame.currentFrame = frameCounter;
            [signalFrame.peakMatrix,signalFrame.powerSpectrumThreshold] = PeakDetection(signalFrame,samplingRate,NFFT,DEBUG);
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

        %if DEBUG == 1
        %    organizedTracks = PlotTracks(frameArray,timeInstants);
        %end
    
    fprintf('\n\n------- SINUSOIDAL ANALYSIS FINISHED ------\n\n');
