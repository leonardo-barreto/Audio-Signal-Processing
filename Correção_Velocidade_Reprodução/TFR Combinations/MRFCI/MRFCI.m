function [TFR_comb, STFT_tensor, time_combined, freq_combined] = MRFCI(varargin) 
% Description:
% 
%     This function computes a combination of different TFRs using directional 
%     information obtained via the structure tensor method. This version 
%     computes the tensor which contains the TFRs to be combined in blocks, 
%     for memory saving.
%
% Input Arguments:
%
%     x - mono input audio signal to be analyzed
% 
%     fs - sampling rate of x (48kHz as default)
% 
%     Nf - array with number of frequency bins in ascending order. The first
%     element is used for the structure tensor
% 
%           ex: Nf = [1024 2048 4096]
% 
%     block_size - analysis block size in seconds (higher block sizes require 
%     more memory)
% 
%     FChT_flags - 1 for a set of Fan-chirp instances and 0 for spectrograms
% 
%           ex: FChT_flags = [0 1 1]
% 
%     asym_flags - 1 for asymmetric windows and 0 to use
%     symmetric ones
% 
%           ex: asym_flags = [0 0 1]
% 
%     Nf_structure_tensor - analysis window length used to compute the
%     spectrogram used in the structure tensor
% 
%     hop - hop used for all TFR, in samples (256 as default).
% 
%     N_alphas - number of alphas used for STFChTs (only alphas >= 0)
% 
%     C_limits - points of reference for the linear combination procedure
% 
%     anisotropy_thr - threshold used for computing anisotropy measure
% 
%     sigma_t - standard deviation in time axis of the two-dimensional 
%     Gaussian filter G
% 
%     sigma_f - standard deviation in frequency axis of the two-dimensional 
%     Gaussian filter G
%
% Outputs:
% 
%     TFR_comb - combined time-frequency representation
%     
%     X_plot - spectrogram computed with analysis window of length plot_ind
%     
%     time_combined - array containing the time stemp of each frame
%     
%     freq_combined - array containing the analog frequency bins
% 
% This function was created by Mauricio do V. M. da Costa 
% (mauricio.costa@smt.ufrj.br).

x = varargin{1};
fs = varargin{2};
Nf = sort(varargin{3});
block_size = varargin{4};
FChT_flags = varargin{5};
asym_flags = varargin{6};
Nf_structure_tensor = varargin{7};
hop = varargin{8};
N_alphas = varargin{9};
C_limits = varargin{10};
range = varargin{11};
sigma_t = varargin{12};
sigma_f = varargin{13};

%% COMPUTING STRUCTURE TENSOR

% Converting G parameters
sigma_m = sigma_t*fs / hop;
sigma_k = sigma_f*Nf_structure_tensor / fs;

Ngauss = 3*max([sigma_m, sigma_k]); % Covers 3 std for each side

% Computing the STFT for the ST
[STFT, freq_struct_tens, time_struct_tens] = spectrogram(x, hanning(Nf_structure_tensor), ...
                                                                          Nf_structure_tensor - hop, Nf_structure_tensor, fs);

% X_structure_tensor = compress(abs(STFT));
% [C, thetas] = structure_tensor(X_structure_tensor, Ngauss, sigma_m, sigma_k, .3);

X_structure_tensor = compress_dB_norm(abs(STFT).^2, range);
[C, v1, v2] = structure_tensor_v2(X_structure_tensor, Ngauss, sigma_m, sigma_k);

C = C/max(max(C));

% %Plotting structure tensor vectors - - - -
% thetas = atan(v2./v1);
% plot_quiver(X_structure_tensor, freq_struct_tens, time_struct_tens, C, thetas);
% return
% %Plotting structure tensor vectors - - - -

% Computing alphas
k_array = 0:length(freq_struct_tens(:))-1;
% alphas = tan(thetas)*fs ./ (hop*repmat(k_array(:), 1, length(time_struct_tens)));
alphas = (v2./v1)*fs ./ (hop*repmat(k_array(:), 1, length(time_struct_tens)));
alphas(1,:) = -10^11;

alphas(isnan(alphas)) = 0;

clear v1 v2 % saving memory space

% Interpolating alphas and C
time_combined = Nf(end)/(2*fs) : hop/fs : ((floor((length(x) - Nf(end)) / hop) + 1)*hop + Nf(end)/2 - 1)/fs;
freq_combined = (0:Nf(end)/2)*fs/Nf(end);

[m_grid_orig, k_grid_orig]  = meshgrid(time_struct_tens, freq_struct_tens);
[m_grid_final, k_grid_final] = meshgrid(time_combined, freq_combined);

alphas_interp = interp2(m_grid_orig, k_grid_orig, alphas, m_grid_final, ...
                                                                    k_grid_final, 'linear', 0);

C_interp = interp2(m_grid_orig, k_grid_orig, C, m_grid_final, k_grid_final, ...
                                                                                        'linear', 0);

% Processing combination in blocks
initial_samples_inds_frames = 1 : hop : ((floor((length(x) - Nf(end)) / hop) + 1)*hop - 1);
final_samples_inds_frames = Nf(end) : hop : ((floor((length(x) - Nf(end)) / hop) + 1)*hop + Nf(end) - 1);

n_blocks = ceil(time_combined(end)/block_size);
n_frames_per_block = floor(length(time_combined)/n_blocks);

TFR_comb = zeros(Nf(end)/2 + 1, length(time_combined));
STFT_tensor = zeros(Nf(end)/2 + 1, length(time_combined), length(Nf));

for block_ind = 1 : n_blocks
%     disp(['Computing block ' num2str(block_ind)])

    if block_ind == 1
        initial_frame = 1;
        initial_sample = 1;
        final_frame = n_frames_per_block*block_ind;
        final_sample = final_samples_inds_frames(final_frame);
    elseif block_ind == n_blocks
        initial_frame = 1 + n_frames_per_block*(block_ind - 1);
        initial_sample = initial_samples_inds_frames(initial_frame);
        final_frame = size(TFR_comb,2);
        final_sample = length(x);
    else
        initial_frame = 1 + n_frames_per_block*(block_ind - 1);
        initial_sample = initial_samples_inds_frames(initial_frame);
        final_frame = n_frames_per_block*block_ind;
        final_sample = final_samples_inds_frames(final_frame);
    end

    x_block = x(initial_sample : final_sample);
    C_block = C_interp(:, initial_frame : final_frame);
    alphas_block = alphas_interp(:, initial_frame : final_frame);

    if sum(abs(x_block)) > 0
        [TFR_comb(:, initial_frame : final_frame), STFT_tensor(:, initial_frame : final_frame, :)] = ...
                                                    lin_comb_mult_res(x_block, fs, Nf, FChT_flags, asym_flags, ...
                                                                C_block, alphas_block, hop, N_alphas, C_limits);
    else
        TFR_comb(:, initial_frame : final_frame)       = 0;
        STFT_tensor(:, initial_frame : final_frame, :) = 0;
    end
end
end

function [TFR_comb, STFT_tensor] = lin_comb_mult_res(varargin)
x = varargin{1};
fs = varargin{2};
Nf = sort(varargin{3});
FChT_flags = varargin{4};
asym_flags = varargin{5};
C = varargin{6};
alphas = varargin{7};
hop = varargin{8};
N_alphas = varargin{9};
C_limits = varargin{10};

%%  PARAMETERS CONFIGURATION * * * * * * * * * * * * * * * * * * * *

% Computing initial alphas
alpha_max = 2*fs/Nf(end);
alphas_ind = 0 : N_alphas;
initial_alphas = tan(alphas_ind * atan(alpha_max)/ N_alphas);
initial_alphas = [-initial_alphas(end:-1:2) initial_alphas];

if size(C_limits,1) == 1
    C_limits = [C_limits ; 1:length(C_limits)];
end

%% COMPUTING TFRs * * * * * * * * * * * * * * * * * * *  

% Computing the STFT with the shortest analysis window - - -
w = hanning(Nf(1));% symetric analysis window
[X_short, Freqs_short, Time_short] = spectrogram(x, w, Nf(1) - hop, Nf(1), fs);
X_short = abs(X_short);

% Interpolating X_short
time_combined = Nf(end)/(2*fs) : hop/fs : ((floor((length(x) - Nf(end)) / hop) + 1)*hop + Nf(end)/2 - 1)/fs;
freq_combined = (0:Nf(end)/2)*fs/Nf(end);

[m_grid_orig, k_grid_orig]  = meshgrid(Time_short, Freqs_short);
[m_grid_final, k_grid_final] = meshgrid(time_combined, freq_combined);

X_short_interp = interp2(m_grid_orig, k_grid_orig, X_short, m_grid_final, ...
                                   k_grid_final, 'linear', 0);

% Initializing tensor for all TFRs to be combined
TFRs_tensor = zeros(length(freq_combined), length(time_combined), ...
                             length(initial_alphas), length(Nf));

STFT_tensor = zeros(length(freq_combined), length(time_combined), length(Nf));

% Compute the TFR for all layers
for ind_Nf = 1:length(Nf)

    if asym_flags(ind_Nf)
        w = asym_hanning(Nf(ind_Nf), Nf(ind_Nf)/2);
    else
        w = hanning(Nf(ind_Nf));
    end
    
	[X, Freq, Time] = spectrogram(x, w, Nf(ind_Nf) - hop, Nf(ind_Nf), fs);
    
    % Interpolation
    [m_grid_orig, k_grid_orig] = meshgrid(Time, Freq);
    X_interp = interp2(m_grid_orig, k_grid_orig, abs(X), m_grid_final, ...
                                k_grid_final, 'linear', 0);                    
                            
    STFT_tensor(:, :, ind_Nf) = X_interp;
%     STFT_tensor(:, :, ind_Nf) = X_interp.^2;

    if FChT_flags(ind_Nf) % Will have fan-chirp instances
        % Compute the instances of STFChTs
        for alphas_ind = 2 : length(initial_alphas) - 1
            if initial_alphas(alphas_ind) ~= 0
                [STFChT, ~, ~] = fanchirpogram(x, w, hop, fs, initial_alphas(alphas_ind));

                STFChT_interp = interp2(m_grid_orig, k_grid_orig, abs(STFChT), ...
                                                    m_grid_final, k_grid_final, 'linear', 0);
                
                TFRs_tensor(:, :, alphas_ind, ind_Nf) = STFChT_interp;
            else
                TFRs_tensor(:, :, alphas_ind, ind_Nf) = X_interp;
            end
        end
        % STFT is used as the first and the last layers
        TFRs_tensor(:, :, 1, ind_Nf) = X_short_interp;
        TFRs_tensor(:, :, end, ind_Nf) = X_short_interp;
    else % Will have only DFT instances
        for alphas_ind = 1 : length(initial_alphas)
            TFRs_tensor(:, :, alphas_ind, ind_Nf) = X_interp;
        end
    end
end

%% ENERGY LEVELLING

energy_ref = sum(sum(X_short_interp.^2));

for ind_Nf = 1 : size(TFRs_tensor, 4)
    for RTF_ind = 1 : size(TFRs_tensor, 3)
        energy = sum(sum(TFRs_tensor(:, :, RTF_ind, ind_Nf).^2));
        TFRs_tensor(:, :, RTF_ind, ind_Nf) = TFRs_tensor(:, :, RTF_ind, ind_Nf)*...
                                                                          sqrt(energy_ref / energy);
    end
end

X_short_interp = X_short_interp.^2;
TFRs_tensor = TFRs_tensor.^2;

%% COMBINING

% Computing the tensor containing the combination weights

TFR_comb = zeros(size(TFRs_tensor,1), size(TFRs_tensor,2));

weights_tensor_alphas = zeros(size(TFRs_tensor, 1), size(TFRs_tensor, 2), size(TFRs_tensor, 3));

for ind_alphas = 1 : size(TFRs_tensor,3) - 1

    aux_alphas_matrix1 = zeros(size(TFR_comb));
    aux_alphas_matrix2 = zeros(size(TFR_comb));    
    
    inds_matrix = (alphas > initial_alphas(ind_alphas)) & ...
                        (alphas <= initial_alphas(ind_alphas + 1));
    
    alpha_values = alphas(inds_matrix);

    
    weight_alpha_2 = (alpha_values - initial_alphas(ind_alphas)) / ...
                            (initial_alphas(ind_alphas + 1) - initial_alphas(ind_alphas));
    
    weight_alpha_1 = 1 - weight_alpha_2;
    
    aux_alphas_matrix1(inds_matrix) = weight_alpha_1;
    aux_alphas_matrix2(inds_matrix) = weight_alpha_2;
    
    weights_tensor_alphas(:, :, ind_alphas) = aux_alphas_matrix1 + ...
                                                      weights_tensor_alphas(:, :, ind_alphas);
    
    weights_tensor_alphas(:, :, ind_alphas + 1) = aux_alphas_matrix2 + ...
                                                weights_tensor_alphas(:, :, ind_alphas + 1);
end

for ind_C = 1 : size(C_limits, 2) - 1
    
    C_matrix1 = zeros(size(TFR_comb));    
    C_matrix2 = zeros(size(TFR_comb));
    
    inds_matrix = (C > C_limits(1, ind_C)) & ...
                        (C <= C_limits(1, ind_C + 1));
    
    C_values = C(inds_matrix);

    weight_C_2 = (C_values - C_limits(1, ind_C)) / (C_limits(1, ind_C + 1) - ...
                                                                            C_limits(1, ind_C));
                                                                            
    weight_C_1 = 1 - weight_C_2;
    
    C_matrix1(inds_matrix) = weight_C_1;
    C_matrix2(inds_matrix) = weight_C_2;

    aux_matrix_1 = C_matrix1 .* sum(weights_tensor_alphas .* ...
                                                        TFRs_tensor(:, :, :, C_limits(2, ind_C)), 3);

    aux_matrix_2 = C_matrix2 .* sum(weights_tensor_alphas .* ...
                                                            TFRs_tensor(:, :, :, C_limits(2, ind_C + 1)), 3);

    TFR_comb(inds_matrix) = aux_matrix_1(inds_matrix) + ...
                                                                aux_matrix_2(inds_matrix);

end

% Filling TFR where C <= C_limits(1)
TFR_comb(C <= C_limits(1, 1)) = X_short_interp(C <= C_limits(1, 1));

% Filling TFR where C > C_limits(end)
aux_matrix = sum(weights_tensor_alphas.*TFRs_tensor(:, :, :, end), 3);
TFR_comb(C > C_limits(1, end)) = aux_matrix(C > C_limits(1, end));

% Filling TFR in the attacks
TFR_comb(alphas < initial_alphas(1)) = X_short_interp(alphas <= ...
                                                                                    initial_alphas(1));
TFR_comb(alphas >= initial_alphas(end)) = X_short_interp(alphas >= ...
                                                                                  initial_alphas(end));
end