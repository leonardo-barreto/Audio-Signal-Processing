% This script performs tests using the MRFCI combination method and plots figures.

% - - - - - - Setting Paths - - - - - - 

if isunix
    figsPath = path_check('./audio/figs/');
else
    figsPath = path_check('.\audio\figs\');
end

%addpath audio
%addpath tests

% - - - Plotting parameters - - -
plot_ind = []; % Chooses the Nf to plot the spectrogram
redLines = []; % Time instants, in s, for ploting vertical-dashed lines

%% - - - - - - Computing MRFCI TFR - - - - - -

% Computing combination procedure in blocks
    [TFR_comb, STFT_tensor, Time_combined, Freq_combined] = MRFCI(x, fs, Nf, block_size,...
        FChT_flags, asym_flags, Nf_structure_tensor, hop, N_alphas, C_limits, range, sigma_t, sigma_f);

%% - - - - - - Plotting individual TFRs - - - - - -

     %dB

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
    %imagesc(Time_combined, Freq_combined, 10*log10(TFR_comb));
    set(gca,'YDir','normal');
    ylim(y_lim_vec);
    colormap(1-gray);
    xlabel('Tempo (s)');
    ylabel('Frequencia (Hz)');
    tit = ['TFRs_' signal_name '_st_lin_comb_Nalphas_' num2str(N_alphas)];

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