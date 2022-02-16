function [freqComponents, timeInstants, TFRepresentation] = TFAnalysis_MRFCI(signal_name)
    
    %   This function makes a time-frequency analysis of a signal, using the MRFCI method.

    if isunix
        addpath ./audio
        dirbar = '/';
    else
        addpath .\audio
        dirbar = '\';
    end
    

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        fs = 44100;
        hop = 128;

        block_size = 2; % Analysis block size in seconds

        % Synthesis TFR - - - - - - - - - - - - - - - - - - -
        N_alphas = 7; % Number of alphas (I, in the paper)

        Nf = [1024 2048 3072 4096]; % Window sizes
        asym_flags = [0 0 0 1]; % Indicate which Nf will be used to compute Fan-chirp instances
        FChT_flags = [0 1 1 1]; % Indicate which Nf will be used to compute Fan-chirp instances

        % Structure tensor - - - - - - - - -

        C_limits = [0 .25 .75 1 ; 1 2 2 3]; % Defines the points of transition for the lambda weights
        range = 60; % range in dB for the ST

        % G parameters
        sigma_t = Nf(end)/(4*fs); % in s
        sigma_f = 100;% in Hz

        Nf_structure_tensor = Nf(1); % To compute the structure tensor

        % Compression
        plot_range = 80; %dB


    %% - - - - - - Input Reading - - - - - - 
        
        [data, fs_orig] = audioread([signal_name]);

        [filepath,signal_name,ext] = fileparts(signal_name);

        if size(data,2) > 1
            x = mean(data.');
        else
            x = data;
        end
        x = x(:);

        if fs ~= fs_orig
            x = resample(x, fs, fs_orig);
        end
        %x = x(1: 15*fs);

        y_lim_vec = [0 5000];


    %% - - - - - - -  TFR computation - - - - - - -

        main_MRFCI;
            
        %Gathering TFRs
        freqComponents = Freq_combined';
        timeInstants = Time_combined;
        
        TFRepresentation = TFR_comb;  
        %TFRepresentation = compress_dB_norm(TFRepresentation, plot_range);
        %TFRepresentation = 10*log10(TFRepresentation);    

    %% - - - - - - -  Plotting - - - - - - -  
        
        %PlotSpectrogram_ylin(freqComponents,timeInstants,[10-plot_range 10],10*log10(TFRepresentation));
end
