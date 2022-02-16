function [freqComponents, timeInstants, TFRepresentation] = TFAnalysis_FLS(signal_name)
    
    %   This function makes a time-frequency analysis of a signal, using the FLS method.

    if isunix
        addpath ./audio
        dirbar = '/';
    else
        addpath .\audio
        dirbar = '\';
    end
    

    %% - - - - - - Parameters configuration - - - - - - 

        % Input Signal params
        resampling_factor = 1; % 1, 1/2 or 1/4
        fs = 44100*resampling_factor;

        % Standard setup
        NFFT = 4096*resampling_factor;
        hop_DFT = 0.0029; % Final time interval between frames in seconds (interp after all computations)
        % hop_DFT = 0.001; % Final time interval between frames in seconds (interp after all computations)

        % % To plot 'sin_1khz'
        % NFFT = 4096*2; % use this to plot 'sin_1khz'

        % % To plot 'impulse'
        % NFFT = 4096/4;
        
        % % Window lengths
        % N_w = [1, 2, 3, 4]*1024*resampling_factor;
        % % N_w = [1024 4096]/4;
        % overlap_short = N_w(1)*0.75; % overlap to be used initially in the spectrograms
        % % overlap_percentage = 0.75; % overlap to be used initially in the spectrograms
        % hop_DFT = 0.005; % Final time interval between frames in seconds (interp after all computations)
        % % hop_DFT = 0.001; % Final time interval between frames in seconds (interp after all computations)
        % % final_time_scale = 0:hop_DFT:length(x)/fs;
        % 
        % % 2D Analysis window
        % time_span = 0.05; % time span of the 2D window in seconds
        % freq_span = 1722; % frequency span of the 2D window in Hz
        % size_W_t_f = [freq_span , time_span]; % In terms of Hz and s

        % Window lengths (in ascendent order!!!)
        N_w = [1, 2, 3, 4]*1024*resampling_factor;
        samp_hop = ceil(hop_DFT*fs);
        %samp_hop = 128;
        overlap_short = N_w(1) - samp_hop;

        % 2D Analysis window
        time_span = 0.03; % time span of the 2D window in seconds
        Ncols_DFT = ceil(time_span/hop_DFT);
        Nrows_DFT = 3*N_w(end)/N_w(1);
        size_W_m_k = [Nrows_DFT Ncols_DFT];

        % size_w_E = [Nrows_DFT/3   2*Ncols_DFT]; %%%%%%%%%%

        % plot_window_E(size_w_E, 'window_E.png', 40)
        % plot_window_S(size_w_S, 'window_S.png', 40)
        % 
        % return

        eta = 20;

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

        sig_len = length(x)/fs;

        y_lim_vec = [0 5000];


    %% - - - - - - -  TFR computation - - - - - - - 
 
        main_FLS;
                
        %Gathering TFRs
        freqComponents = f;
        timeInstants = t;
        
        TFRepresentation = combined_TFR_HLS;   
        %TFRepresentation = compress_dB_norm(TFRepresentation, plot_range);
        %TFRepresentation = 10*log10(TFRepresentation);    

    %% - - - - - - -  Plotting - - - - - - -  
        
        %PlotSpectrogram_ylin(freqComponents,timeInstants,[10-plot_range 10],10*log10(TFRepresentation));

end
