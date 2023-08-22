clear;
%% - - - - - - - - - - - - - - - INPUT PARAMETERS - - - - - - - - - - - - - - - 
    %
    % This section is where you should edit your desired parameters.
    % It is not recommended to edit anything further than this section, unless you run into unexpected errors.
    %
    % Parameters overview:
    % - Audio file name
    % - Fixed sampling rate (if pre-resampling is needed)
    % - Time-frequency analysis methods
    % - Time-varying resampling parameters

    fileName = 'paulistana3_5s.wav';
    %fileName = 'source_Midsum.wav';
    %fileName = 'Mix1.wav';

    fs = 44100;     % Desired fixed sampling rate (if pre-resampling the original signal is needed)

    % TFR Method
    method_name = {'STFT', 'CQT', 'SWGM', 'FLS', 'FEMD', 'MRFCI'};      % TFR Methods available
    method_flags = [1 0 1 0 0 1];                                       % Flags for which method(s) will be enabled
    backwardsFlag = 0;                                                  % Flag for backwards sinusoidal tracking
    methods_enabled = find(method_flags);                

    % HPSS Options
    HPSS_on = [];                   % Harmonic-percussive separation flag 
                                        % For turning on: 'HPSS'. 
                                        % For turning off: [] (empty variable).
    if HPSS_on
        nFilterSS = 71;             % SS filter filter size: must be odd
        nFilterTr = 71;             % Transient filter size: must be odd
        nIter = 3;                  % No. of HPSS iterations
        HPSS_method = 'median';     % 'median' or 'SSE'
        kernel_option = 'normal';   % 'normal' or 'relaxed'
        HPSS_options = struct('nFilterSS',nFilterSS,'nFilterTr',nFilterTr,'nIter',nIter,'method',HPSS_method,'kernel_option',kernel_option);
    else
        HPSS_options = [];
    end

    % Resampling parameters
    filterCoeffs = 1024;    % Number of coefficients for the resampling sinc filter
    pitchOffset = 0;        % Global pitch offset (PERCENTAGE - leave it 0 for no pitch offset at all)

    % Flags and general
    plot_enable = 1;
    print_figures = 1;
    energy_ref_method = 1;  % index of method that will be used as reference energy for plots (guide in method_name)
    plot_range = 100;       % dB - Power range for plotting
    %plot_max = 10;          % dB - Max plotting power

%% - - - - - - Reading input audio file and setting path for figures - - - - - -

    [inputSignal,signalName] = ReadAudioFile(fileName, fs);

    if isunix
        addpath ./audio_src
        dirbar = '/';
        figsPath = path_check('./figures_out/SpdCorr_Experiment/');
        figsPath = path_check([figsPath '/' signalName]);
    else
        addpath .\audio_src
        dirbar = '\';
        figsPath = path_check('.\figures_out\SpdCorr_Experiment\');
        figsPath = path_check([figsPath '\' signalName]);
    end
    

%% - - - - - - Sinusoidal Analysis - - - - - -  

    fprintf('\n------- SINUSOIDAL ANALYSIS ------\n');
    TFR = cell(length(methods_enabled,1));
    signalTracks = cell(length(methods_enabled,1));
    TFParams = cell(length(methods_enabled,1));
    for i = methods_enabled
        if HPSS_on
            [TFR{i},TFParams{i},signalTracks{i}] = SinusoidalAnalysis(inputSignal,fs,method_name{i},backwardsFlag,HPSS_on,HPSS_options);
        else
            [TFR{i},TFParams{i},signalTracks{i}] = SinusoidalAnalysis(inputSignal,fs,method_name{i},backwardsFlag);
        end
    end

    % Making first frame zero
    for i = methods_enabled
        TFParams{i}.timeInstants = TFParams{i}.timeInstants-TFParams{i}.timeInstants(1);
    end

    % Energy normalization
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

    % Plotting TFRs and tracks
    if plot_enable
        for i = methods_enabled
            plot_method_name = method_name{i};
            if HPSS_on
                plot_method_name = [method_name{i} '-SS'];
            end
            plot_max = max(max(10*log10(TFR{i})));
            PlotSpectrogram_ylin(TFParams{i}.freqComponents,TFParams{i}.timeInstants,[plot_max-plot_range plot_max],10*log10(TFR{i}));
            title(sprintf('RTF (%s)',plot_method_name));

            if print_figures
                tit = [signal_name '_RTF_' plot_method_name];
                tit(tit=='.') = '_'; tit(tit==' ') = '';
                figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
                figFileName = [figsPath tit];
                formatFig(gcf, figFileName, 'en', figProp);
                close;
            end

            organizedTracks = PlotTracks(signalTracks{i},TFParams{i}.timeInstants);
            title(sprintf('Trilhas senoidais (%s)',plot_method_name));

            if print_figures
                tit = [signal_name '_tracks_' plot_method_name];
                tit(tit=='.') = '_'; tit(tit==' ') = '';
                figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
                figFileName = [figsPath tit];
                formatFig(gcf, figFileName, 'en', figProp);
                close;
            end

        end
    end

%% - - - - - - Speed Correction - - - - - -  

    resamplingFactors = cell(length(methods_enabled,1));
    outputSignal = cell(length(methods_enabled,1));
    for i = methods_enabled
        fprintf('\n------- PITCH VARIATION CURVE EXTRACTION ------\n');
        fprintf('\nPVC extraction in progress...\n');
        resamplingFactors{i} = ExtractPVC(signalTracks{i},TFParams{i}.totalFrames);
        resamplingFactors = resamplingFactors+(pitchOffset/100);
        fprintf('PVC extraction done.\n');
        fprintf('\n------- TIME-VARYING RESAMPLING ------\n');
        fprintf('\nTime-varying resample in progress (this might take some time)...\n');
        outputSignal{i} = TimeVarying_Resample(inputSignal,TFParams{i},resamplingFactors{i},filterCoeffs);
        fprintf('Time-varying resample done.\n');
    end

    % Plotting PVC
    if plot_enable
        for i = methods_enabled
            plot_method_name = method_name{i};
            if HPSS_on
                plot_method_name = [method_name{i} '-SS'];
            end
            PlotPVC(resamplingFactors{i},TFParams{i}.timeInstants);
            title(sprintf('Curva de desvio de pitch (%s)',plot_method_name));

            if print_figures
                tit = [signal_name '_PVC_' plot_method_name];
                tit(tit=='.') = '_'; tit(tit==' ') = '';
                figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
                figFileName = [figsPath tit];
                formatFig(gcf, figFileName, 'en', figProp);
                close;
            end
        end
    end

load handel
x = y;
load laughter
lolSize = min(length(y),length(x));
sound(y(1:lolSize)/8+x(1:lolSize)/10,Fs)
clearvars x y Fs
