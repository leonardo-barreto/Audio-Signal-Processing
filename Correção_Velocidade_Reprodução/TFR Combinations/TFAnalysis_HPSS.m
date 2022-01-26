function [TFParams] = TFAnalysis_HPSS(signal_name)
    
    %   This function makes a time-frequency analysis of a signal, using auxiliary functions for modularity.
    %
    %   1st - TFR
    %

    DEBUG = 1;

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        fs = 44100;
        hop = 1024;

        block_size = 2; % Analysis block size in seconds

        % Synthesis TFR - - - - - - - - - - - - - - - - - - -
        N_alphas = 7; % Number of alphas (I, in the paper)

        Nf = [2048 4096 8192]; % Window sizes
        asym_flags = [0 0 0]; % Indicate which Nf will be used to compute Fan-chirp instances
        FChT_flags = [0 0 1]; % Indicate which Nf will be used to compute Fan-chirp instances

        % Structure tensor - - - - - - - - -

        C_limits = [0 .25 .75 1 ; 1 2 2 3]; % Defines the points of transition for the lambda weights
        range = 60; % range in dB for the ST

        % G parameters
        sigma_t = Nf(end)/(4*fs); % in s
        sigma_f = 100;% in Hz

        Nf_structure_tensor = Nf(3); % To compute the structure tensor

        % Compression
        plot_range = 80; %dB



    %% - - - - - - Input Reading - - - - - - 
    fprintf('\n\n------- READING INPUT ------\n\n')
        
        [data, fs_orig] = audioread([signal_name]);

        if size(data,2) > 1
            x = mean(data.');
        else
            x = data;
        end
        x = x(:);

        %x = resample(x, fs, fs_orig);
        %x = x(1: 15*fs);
        inputSignal = x;

        ylim_vec = [0 4000];


    %% -------------------------------------- TFR stage -------------------------------------------
    fprintf('\n\n------- TIME-FREQUENCY ANALYSIS STARTED ------\n\n');
        % MRFCI
            main_MRFCI_HPSS;
                
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
        

        %Outputting general parameters.
        TFParams = {};
        TFParams.signalSize = length(inputSignal);
        TFParams.samplingRate = samplingRate;
        TFParams.timeInstants = timeInstants;
        TFParams.fftPoints = fftPoints/2;
        TFParams.totalFrames = totalFrames;
        TFParams.hopSize = hopSize;

        
        fprintf(' Total frames: %i\n',TFParams.totalFrames);
        fprintf(' Number of frequency bins: %i\n',TFParams.fftPoints);
        fprintf('\nTime-Frequency Representation Finished.\n');

        %if DEBUG == 1 %call plot
            %PlotSpectrogram(freqComponents,timeInstants,powerMatrixDB);
        %end

    
fprintf('\n\n------- TIME-FREQUENCY ANALYSIS FINISHED ------\n\n');
