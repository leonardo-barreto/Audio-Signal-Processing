clear;
%   This function makes a time-frequency analysis of a signal.

if isunix
    addpath ./audio
    dirbar = '/';
else
    addpath .\audio
    dirbar = '\';
end

%% - - - - - - Input Parameters - - - - - - 
    
    signal_name = 'violin_drum2.wav';
    fs = 44100;

    % Method
    method_name = {'STFT', 'MRFCI', 'FLS'}; % TFR Methods available
    method_flags = [1 0 0]; % Which method will be enabled

    % Plotting parameters
    plot_enable = 1; % 1 enables plotting, 0 disables

    energy_ref_method = 2; % index of method that will be used as reference energy for plots (guide in method_name)
    plot_range = 80; % dB - Power range for plotting
    plot_max = 10; % dB - Max plotting power

%% - - - - - - Input Reading - - - - - -  
    
    [data, fs_orig] = audioread([signal_name]);

    %[filepath,signal_name,ext] = fileparts(signal_name);

    if size(data,2) > 1
        x = mean(data.');
    else
        x = data;
    end
    x = x(:);

    if fs ~= fs_orig
        x = resample(x, fs, fs_orig);
    end 

%% - - - - - - -  TFR computation - - - - - - - 

    fprintf('\nSignal under analysis: %s\n\n', signal_name)

    TFR_function = {@TFAnalysis_STFT @TFAnalysis_MRFCI @TFAnalysis_FLS};
    methods_enabled = find(method_flags);

    f = {};
    t = {};
    TFR = {};

    for i = methods_enabled
        fprintf('---- %s TFR ----\n',method_name{i});

        [f{i},t{i},TFR{i}] = TFR_function{i}(x, fs);

        fprintf('   Total frames: %i\n',length(t{i}));
        fprintf('   Frequency bins (N_w/2 + 1): %i\n\n',length(f{i}));
    end

%% - - - - - - - Energy normalization - - - - - - -

    if (numel(methods_enabled) > 1)
        if (~isempty(find(energy_ref_method == methods_enabled)))
            ref_energy = sum(sum(TFR{energy_ref_method}));
            for i = methods_enabled
                if i ~= energy_ref_method
                    TFR{i} = ref_energy.* TFR{i}/sum(sum(TFR{i}));
                end
            end
        else
            error('Invalid method for energy reference. Choose one of the active options.')
        end
    end


%% - - - - - - - Plotting - - - - - - -
    if plot_enable
        for i = methods_enabled
            PlotSpectrogram_ylin(f{i},t{i},[10-plot_range 10],10*log10(TFR{i}));
            title(sprintf('Espectrograma %s', method_name{i}))
            %set(gca, 'FontSize', 30, 'yscale', 'log')
            %ylim([0 12000])
        end
    end


%% ------ TEMPORARY VARIABLE CLEAR ------
    clearvars -except f t TFR 
    
