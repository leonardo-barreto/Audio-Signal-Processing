%TFAnalysis_all;
clear;
%clearvars -except signals_names signal_name;

if isunix
    addpath ./audio
    dirbar = '/';
    figsPath = path_check('./Temporary_Figures/spd_correction/');
else
    addpath .\audio
    dirbar = '\';
    
end

%% - - - - - - Input Parameters - - - - - - 
    
    %signal_name = 'source_Midsum.wav';
    signal_name = 'paulistana3_5s.wav';
    %signal_name = 'Mix1.wav';

    fs = 44100;

    % TFR Method
    method_name = {'STFT', 'FLS', 'MRFCI'}; % TFR Methods available
    method_flags = [1 1 1]; % Which method will be enabled
    methods_enabled = find(method_flags);

    plot_enable = 1;
    print_figures = 1;
    plot_range = 80;

    energy_ref_method = 3; % index of method that will be used as reference energy for plots (guide in method_name)
    %plot_range = 100; % dB - Power range for plotting
    %plot_max = 10; % dB - Max plotting power

%% - - - - - - Input Reading - - - - - -  
    
    [data, fs_orig] = audioread([signal_name]);

    [filepath,signal_name,ext] = fileparts(signal_name);
    signal_name = char(signal_name);

    if size(data,2) > 1
        x = mean(data.');
    else
        x = data;
    end
    x = x(:);

    if fs ~= fs_orig
        x = resample(x, fs, fs_orig);
    end 


TFR = {};
signalTracks = {};
TFParams = {};


for i = methods_enabled
    [TFR{i},signalTracks{i},TFParams{i}] = SinusoidalAnalysis(x,fs,method_name{i});
end

%% - - - - - - Making first frame zero
for i = methods_enabled
    TFParams{i}.timeInstants = TFParams{i}.timeInstants-TFParams{i}.timeInstants(1);
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

if plot_enable
    for i = methods_enabled
        plot_max = max(max(10*log10(TFR{i})));
        PlotSpectrogram_ylin(TFParams{i}.freqComponents,TFParams{i}.timeInstants,[plot_max-plot_range plot_max],10*log10(TFR{i}));
        title(sprintf('RTF (%s)',method_name{i}));

        if print_figures
            figsPath = path_check(['.\Temporary_Figures\spd_correction\' signal_name '\']);
            tit = [signal_name '_RTF_' method_name{i}];
            tit(tit=='.') = '_'; tit(tit==' ') = '';
            figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
            figFileName = [figsPath tit];
            formatFig(gcf, figFileName, 'en', figProp);
            close;
        end

        organizedTracks = PlotTracks(signalTracks{i},TFParams{i}.timeInstants);
        title(sprintf('Trilhas senoidais (%s)',method_name{i}));

        if print_figures
            figsPath = path_check(['.\Temporary_Figures\spd_correction\' signal_name '\']);
            tit = [signal_name '_tracks_' method_name{i}];
            tit(tit=='.') = '_'; tit(tit==' ') = '';
            figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
            figFileName = [figsPath tit];
            formatFig(gcf, figFileName, 'en', figProp);
            close;
        end

    end
end

%% ------ TEMPORARY VARIABLE CLEAR ------
%clearvars data dirbar ref_energy energy_ref_method i method_flags method_name methods_enabled plot_enable plot_max plot_range signal_name TFR_function method
%clearvars ans figFileName figProp figsPath nFilterSS nFilterTr nIter print_figures tit
    