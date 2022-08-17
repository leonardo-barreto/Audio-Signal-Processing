%clear variables; clc; close all;

%%%%%%%%%%%%%
%       Setting Paths
%%%%%%%%%%%%%

TFR_files_dir = path_check('/Volumes/HD_MAURICIO/datasets/MedleyDB/MDB-melody-synth/MedleyDB_MRFCI_TFRs_config1');
audio_files_dir = '/Volumes/HD_MAURICIO/datasets/MedleyDB/MDB-melody-synth/audio_melody';

%%%%%%%%%%%%%%%%%%
%       Parameters configuration     %
%%%%%%%%%%%%%%%%%%

read_mdb_synth_config1

%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
%%%%%%%%%%%%%%%%
%       Reading input signal      %
%%%%%%%%%%%%%%%%

wav_file_name = files_in_path(audio_files_dir,{'wav'});

for file_ind = 1:length(wav_file_name)    
    %Load wav file - - - - - - - - - - - - - - - - - - - - -
    [data, fs_orig] = audioread(wav_file_name(file_ind).name);
    
    signal_name = wav_file_name(file_ind).clean_file_name;
    
    disp([signal_name ' file ' num2str(file_ind) ' of ' num2str(length(wav_file_name))])
        
    % Make it mono
    if size(data,2) > 1
        x = mean(data.');
    else
        x = data;
    end
    x = x(:);
    
    % Resampling
    if fs_orig ~= fs
        x = resample(x, fs, fs_orig);
    end    

    % Normalizing x
    x = x/sqrt(mean(x.^2));

    %% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
    %%%%%%%%%%%%%%
    %       Structure Tensor     %
    %%%%%%%%%%%%%%

    % Computing the STFT:------------------------------------------------------
    [STFT, F, time] = spectrogram(x, hanning(Nf),Nf-hop,Nf,fs); 

    X_original = abs(STFT); 

    X = X_original;

    % Computing the structure tensor:------------------------------------------
    X_structure_tensor = compress_dB_norm(X.^2, range);

    % Using imgradientxy:
    [Xn, Xw] = imgradientxy(X_structure_tensor);

    % Computing T11, T21, T12 and 22:
    G = gauss2D(Ngauss, [sigma_n sigma_k]);
    T11 = conv2(Xn.*Xn, G, 'same');
    T21 = conv2(Xw.*Xn, G, 'same'); T12 = T21;
    T22 = conv2(Xw.*Xw, G, 'same');

    % Computing angles and anisotropy measure:-----------------------------
    alphas = zeros(size(X)); C = zeros(size(X));
    for m = 1:size(X,2) %time
        for k = 1:size(X,1) %frequency
            T = [T11(k,m) T12(k,m); T21(k,m) T22(k,m)];
            [eig_vec, eig_val] = eig(T);

            % Angle of orientation:
            alphas(k,m) = atan(eig_vec(2,1)/eig_vec(1,1));

            % Anisotropy measure:
            if X_structure_tensor(k,m) > 0
                C(k,m) = ((eig_val(4) - eig_val(1)) / (eig_val(1) + eig_val(4)))^2;
            end            
        end
    end

    R = fs^2/(hop*(size(X,1)-1)*2)*tan(alphas);

    alphasFChT = R./repmat(F(:),1,size(X,2));
    alphasFChT(1,:) = -10^10;
    alphasEST  = zeros(length(time), NumOfSources); 

    % Gaussian Kernel:
    kernel = @(x)(exp(-.5*x.^2)/sqrt(2*pi));

    for m = 1:size(X,2)
        % Defining PDF:
        if sum(C(:,m))~=0
            pdfest = @(x)((C(:,m))'*kernel((x-alphasFChT(:,m))/h)/h)/sum(C(:,m)); % Weighting by C 
        else
            pdfest = @(x)mean(kernel((x-alphasFChT(:,m))/h)/h);
        end

        [xpdf,ypdf] = fplot(pdfest, [-10,10], 'k', 100);

        % Estimating alpha
        if m == 1
            alphasEST(m,:) = zeros(1,NumOfSources);
        else
            alphasEST(m,:) = find_alphas(xpdf,ypdf, percent_trsh, NumOfSources);
        end        
    end

    alphasEST_orig = alphasEST;
    if filter_alphas_flag
        alphasEST = filter_alphas_median(alphasEST, ord_median);
        alphasEST = filter_alphas(alphasEST, ord);
    end

    alphasEST(abs(alphasEST)<10^-8) = 0;

    % Compute STFChTs using alphas: - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    [X_2, F_2, time_2] = spectrogram(x, hanning(Nf_2), Nf_2 - hop_2, Nf_2, fs); 
    F_final = F_2;
    X_2 = abs(X_2);

    Num_frames_2 = length(time_2); 
    STFChT = zeros(Nf_2/2+1, Num_frames_2, NumOfSources);
    alphasEST_interp = zeros(length(time_2), NumOfSources);

    for alphas_ind = 1:NumOfSources
        % Interpolation of the estimated alphas   
        alphasEST_interp(:,alphas_ind) = interp1(time, alphasEST(:, alphas_ind), time_2, 'pchip');    
        for ii = 1:Num_frames_2
            t0 = (ii-1)*hop_2+1;
            xn = x(t0:t0+Nf_2-1);
            [~,mag_specs] = calculaFChT(xn, fs, Nf_2, 1, alphasEST_interp(ii, alphas_ind));
            STFChT(:,ii,alphas_ind) = abs(mag_specs(1:Nf_2/2+1));
        end
    end

    if combine_original_TFR
        [m_grid, k_grid] = meshgrid(time, F);
        [m_grid_2, k_grid_2] = meshgrid(time_2, F_2);
        X_1 = interp2(m_grid, k_grid, X_original, m_grid_2, k_grid_2, 'linear', 0);
        cat_spectrograms = cat(3, X_1, STFChT);
    else
        cat_spectrograms = STFChT;
    end

    if combine_large_STFT
        cat_spectrograms = cat(3, X_2, cat_spectrograms);
    end

    if combine_extra_large_STFT
        [X_3, F_3, time_3] = spectrogram(x, hanning(Nf_3), Nf_3 - hop_2, Nf_3, fs); 
        X_3 = abs(X_3);

        [m_grid, k_grid] = meshgrid(time_3, F_3);
        [m_grid_2, k_grid_2] = meshgrid(time_2, F_2);
        cat_spectrograms_interp = X_3; clear X_3

        for ind = 1:size(cat_spectrograms,3)
            cat_spectrograms_interp = cat(3, cat_spectrograms_interp, interp2(m_grid_2, k_grid_2, ...
                                                        cat_spectrograms(:,:,ind), m_grid, k_grid, 'linear', 0));
        end
        F_final = F_3;
        cat_spectrograms = cat_spectrograms_interp;
    end

    Time_combined = time_2;
    Freq_combined = F_final;
    
    switch comb_method                
        case 1 % LS
            TFR_FEMD = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
            % Saving TFR in a .mat file
            save([TFR_files_dir '/FEMD_LS/' signal_name '.mat'], 'TFR_FEMD', 'Time_combined', 'Freq_combined');                

        case 2 % SLS
            TFR_FEMD_SLS = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
            % Saving TFR in a .mat file
            save([TFR_files_dir '/FEMD_SLS/' signal_name '.mat'], 'TFR_FEMD_SLS', 'Time_combined', 'Freq_combined');
            
        case 3 % All methods
            % LS - - - - - - - - - 
            zeta = 200; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
            eta = 1; % controls the energy weighting process
            TFR_FEMD_LS = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);

            % Saving TFR in a .mat file
            save([path_check([TFR_files_dir '/FEMD_LS/']) signal_name '.mat'], 'TFR_FEMD', 'Time_combined', 'Freq_combined');

            % SLS - - - - - - - - - 
            zeta = 70; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
            eta = 1; % controls the energy weighting process
            TFR_FEMD_SLS = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);

            % Saving TFR in a .mat file
            save([path_check([TFR_files_dir '/FEMD_SLS/']) signal_name '.mat'], 'TFR_FEMD', 'Time_combined', 'Freq_combined');
    end

end