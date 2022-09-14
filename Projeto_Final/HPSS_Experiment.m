%TFAnalysis_all;

if isunix
    figsPath = path_check('./figures_out/Mix1/HPSS/');
else
    figsPath = path_check('.\figures_out\Mix1\HPSS\');
end

% Method
method_name = {'STFT', 'MRFCI', 'FLS'}; % TFR Methods available
method_flags = [1 1 1]; % Which method will be enabled
methods_enabled = find(method_flags);

nFilterSS = 71; % Must be odd
nFilterTr = 71; % Must be odd

nIter = 5;

method = 'median'; % 'median' or 'SSE'

plot_enable = 0;
print_figures = 0;
plot_range = 100;

TFR_SS = {};
TFR_Tr = {};
TFR_Res = {};

for i = methods_enabled
    [TFR_SS{i}, TFR_Tr{i}, TFR_Res{i}] = Iterative_HPR_Separation(TFR{i}, nFilterSS, nFilterTr, nIter, method);
end

%% - - - - - - Making first frame zero
for i = methods_enabled
    t{i} = t{i}-t{i}(1);
end

if plot_enable
    plot_max = max(max(10*log10(TFR{methods_enabled(1)})));
    for i = methods_enabled

        PlotSpectrogram_ylin(f{i},t{i},[plot_max-plot_range plot_max],10*log10(TFR_SS{i}));
        title(sprintf('Estado Permanente (filtro %ix%i, %s)', nFilterSS,nFilterTr,method_name{i}));

        if print_figures
            tit = ['i_' num2str(nIter) '_' num2str(nFilterSS) '_' num2str(nFilterTr) '_' method_name{i} '_SS_enhanced_'];
            tit(tit=='.') = '_'; tit(tit==' ') = '';
            figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
            figFileName = [figsPath tit];
            formatFig(gcf, figFileName, 'en', figProp);
            close;
        end


        PlotSpectrogram_ylin(f{i},t{i},[plot_max-plot_range plot_max],10*log10(TFR_Tr{i}));
        title(sprintf('Transitorio (filtro %ix%i, %s)', nFilterSS,nFilterTr,method_name{i}));

        if print_figures
            tit = ['i_' num2str(nIter) '_' num2str(nFilterSS) '_' num2str(nFilterTr) '_' method_name{i} '_Tr_enhanced_'];
            tit(tit=='.') = '_'; tit(tit==' ') = '';
            figProp = struct('size', 15,'font','Helvetica','lineWidth',2,'figDim',[1 1 560 420]); % Thesis
            figFileName = [figsPath tit];
            formatFig(gcf, figFileName, 'en', figProp);
            close;
        end

    end
end

%% ------ TEMPORARY VARIABLE CLEAR ------
clearvars data dirbar ref_energy energy_ref_method i method_flags method_name methods_enabled plot_enable plot_max plot_range signal_name TFR_function method
clearvars ans figFileName figProp figsPath nFilterSS nFilterTr nIter print_figures tit
    