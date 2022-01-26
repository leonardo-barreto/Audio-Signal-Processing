% This script performs tests using the MRFCI combination method and plots figures.

% - - - - - - Setting Paths - - - - - - 

if isunix
    figsPath = path_check('./figs/');
else
    figsPath = path_check('.\figs\');
end

addpath audio
%addpath tests

% - - - - - - Parameters configuration - - - - - - 

% Analysis TFR - - - - - - - - - - - - - - - - - - - -

fs = 44100;
hop = 1024;

block_size = 2; % Analysis block size in seconds

% Synthesis TFR - - - - - - - - - - - - - - - - - - -
N_alphas = 7; % Number of alphas (I, in the paper)

% Nf = [1024 2048]; % Window sizes
% asym_flags = [0 0]; % Indicate which Nf will be used to compute Fan-chirp instances
% FChT_flags = [1 1]; % Indicate which Nf will be used to compute Fan-chirp instances

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

% - - - Plotting parameters - - -
plot_ind = []; % Chooses the Nf to plot the spectrogram
y_lim_vec = [0 5000];
redLines = []; % Time instants, in s, for ploting vertical-dashed lines

%% - - - - - - Reading input signal - - - - - -
    % read_input
    [data, fs_orig] = audioread([signal_name '.wav']);

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

    % ------ ORIGINAL READING SCRIPTS FROM MAURICIO ------
        % read_synthetic_signal % proof of concept
        % read_harmonic_sin; %<<<<<
        % read_harmonic_sin_var; %<<<<<
        % compute_chirp; %<<<<<<
        % read_violins_vibrato;
        % read_violins_vibrato_drums2;
        % read_pop; % <<<<<
        % read_violin_vocal; % <<<<<
        % read_opera_vocals;
        % read_operaTenor;
        % read_input;
        % read_hadley;
        % read_violins_piano; % xxxxx
        % read_signals; % xxxxx
        % read_harmonic_sin_1;
        % read_carnaval_percussao
        % % read_piano
        % read_apesar_de_voce
        % read_drum_and_bass
        % read_harmonic_chirp_signal
    

    %% - - - - - - Computing MRFCI TFR - - - - - -

    % Computing combination procedure in blocks
    [TFR_comb, STFT_tensor, Time_combined, Freq_combined] = MRFCI(x, fs, Nf, block_size,...
        FChT_flags, asym_flags, Nf_structure_tensor, hop, N_alphas, C_limits, range, sigma_t, sigma_f);

%% - - - - - - Plotting individual TFRs - - - - - -

    plot_range = 80; %dB

    for ind = plot_ind
        figure;
        imagesc(Time_combined, Freq_combined, compress_dB_norm(STFT_tensor(:,:,ind), plot_range));
        set(gca,'YDir','normal');    
        ylim(y_lim_vec); 
        colormap(1-gray);
        xlabel('Time [s]');
        ylabel('Frequency [Hz]');     

        % To show onset and offset
        if ~isempty(redLines)
            hold on;
            for ii = 1:length(redLines)
                y = 0:100:y_lim_vec(2); 
                x = redLines(ii)*ones(size(y)); 
                plot(x,y,'--r');
            end
            hold off;
        end
        
        tit = ['TFRs_' signal_name '_spectrogram_' num2str(Nf(ind))];

        tit(tit=='.') = '_'; tit(tit==' ') = '';

        figProp = struct('size', 35,'font','Times','lineWidth',2,'figDim',[1 1 900 500]); % Thesis
        figFileName = [figsPath tit];
        formatFig(gcf, figFileName, 'en', figProp);
    end

    figure;
    imagesc(Time_combined, Freq_combined, compress_dB_norm(TFR_comb, plot_range));
    set(gca,'YDir','normal');
    ylim(y_lim_vec);
    colormap(1-gray);
    xlabel('Tempo (s)');
    ylabel('Frequencia (Hz)');
    tit = ['TFRs_' signal_name '_st_lin_comb_Na_' num2str(N_alphas)];

    % To show onset and offset
    if ~isempty(redLines)
        hold on;
        for ii = 1:length(redLines)
            y = 0:100:y_lim_vec(2); 
            x = redLines(ii)*ones(size(y)); 
            plot(x,y,'--r');
        end
        hold off;
    end
    
    tit(tit=='.') = '_'; tit(tit==' ') = '';

figProp = struct('size', 35, 'font', 'Times', 'lineWidth', 2 , 'figDim', [1 1 900 500]); % Thesis
    figFileName = [figsPath tit];

    formatFig(gcf, figFileName, 'en', figProp);