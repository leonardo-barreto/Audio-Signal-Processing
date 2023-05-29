%   This script makes a time-frequency analysis of a signal.

if isunix
    figsPath = path_check('./figures_out/TFRs/');
    addpath ./audio_src
    dirbar = '/';
else
    figsPath = path_check('.\figures_out\TFRs\');
    addpath .\audio_src
    dirbar = '\';
end

%% - - - - - - Input Parameters - - - - - - 
    
    %signal_name = 'violin_drum2.wav';
    %signal_name = 'MusicDelta_Hendrix_STEM_04.RESYN.wav';
    %signal_name = 'Paulistana1_5s.wav';
    signal_name = 'Mix1.wav';
    %signal_name = 'paulistana3_5s.wav';
    %signal_name = 'MusicDelta_Pachelbel_STEM_04_Cut.RESYN.wav';

    fs = 44100;

    % Method
    method_name = {'STFT', 'MRFCI', 'FLS'}; % TFR Methods available
    method_flags = [0 1 0]; % Which method will be enabled

    % Plotting parameters
    plot_enable = 1; % 1 enables plotting, 0 disables
    print_figures = 1;

    energy_ref_method = 1; % index of method that will be used as reference energy for plots (guide in method_name)
    plot_range = 100; % dB - Power range for plotting
    %plot_max = 10; % dB - Max plotting power

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

%% - - - - - - -  TFR computation - - - - - - - 

    fprintf('\nSignal under analysis: %s\n\n', signal_name)

    TFR_function = {@TFAnalysis_STFT @TFAnalysis_MRFCI @TFAnalysis_FLS};
    methods_enabled = find(method_flags);

    f = {};
    t = {};
    TFR = {};

    for i = methods_enabled
        fprintf('---- %s TFR ----\n',method_name{i});

        [TFR{i},f{i},t{i}] = TFR_function{i}(x, fs);

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
        plot_max = max(max(10*log10(TFR{methods_enabled(1)})));
        for i = methods_enabled

            PlotSpectrogram_ylin(f{i},t{i},[plot_max-plot_range plot_max],10*log10(TFR{i}));
            title(sprintf('Espectrograma %s', method_name{i}))

            %set(gca, 'FontSize', 30, 'yscale', 'log')
            %ylim([0 12000])
            if print_figures
                tit = [signal_name '_RTF_' method_name{i}];
                tit(tit=='.') = '_'; tit(tit==' ') = '';
                figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
                figFileName = [figsPath tit];
                formatFig(gcf, figFileName, 'en', figProp);
                close;
            end
        end
    end


%% ------ TEMPORARY VARIABLE CLEAR ------
%clearvars data dirbar energy_ref_method i method_flags method_name methods_enabled plot_enable plot_max plot_range signal_name
    