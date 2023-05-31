function [STFChT_comb, F_final, time_2] = TFAnalysis_FEMD(x,fs)
    
    %   This function makes a time-frequency analysis of a signal, using the FEMD method.

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -
        fs = 44100;
        Nf = 1024; hop = Nf/4;
        comb_method = 3; % 1- SWGM, 2 - LS, 3 - SLS, 4 - All methods
        combine_original_TFR = 1; % Includes the initial TFR used by the ST for the alphas estimation <<<<<<<<<<<<<
        combine_large_STFT = 1; % Includes the STFT computed with the same window size as the STFChT
        combine_extra_large_STFT = 0; % Includes an STFT computed with a larger window size than the STFChT
        Nf_3 = Nf*8;

        % use_init_comb = 0; % For combining FChT instances before the ST procedure
        % initial_alphas = [-3 3]; % For combining FChT instances before the ST procedure

        % Final TFR - - - - - - - - - - - - - - - - - - -
        resMult = 2; % Multiplies the resolution of the Analysis TFR
        Nf_2 = Nf*resMult;
        hop_2 = hop;

        useTheoAlphas = 0;  % For testing with controlled signals
        % plotTheoAlphas = 0; % For testing with controlled signals
        renyi_ord = 3;% For evaluating the sparsity of final results
        Num_it = 1;   % Number of iterations
        filter_alphas_flag = 1; % Alphas filtering
        ord = 5;% Alphas filtering
        ord_median = 5;% Alphas filtering

        % - - - Structure Tensor - - -
        sigma_n = 1.5; sigma_k = 1.5; Ngauss = 9; % G parameters
        h = 0.1; % Standard deviation for the distribution
        percent_trsh = 0.05; % Percent of the highest peak to consider as relevant (GMM)
        epsilon = 1;

        % - - - Local Sparsity Combination - - -
        switch comb_method
            case 2 % LS
                zeta = 202; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
                eta = 1; % controls the energy weighting process
            case 3 % SLS
                zeta = 70; % controls the sparsity weighting process (zeta >= 200 for max(LocalSparsity))
                eta = 1; % controls the energy weighting process
        end

        time_span = 0.03; % time span of the 2D window in seconds
        Ncols_DFT = ceil(fs*time_span/hop);
        size_w_S = [21 Ncols_DFT];
        size_w_E = [9  2*Ncols_DFT];

        % - - - SWGM Combination - - -
        beta = .5;

        % EXTRA
        ylim_vec = [0 fs/2]; % in kHz
        NumOfSources = 1; %Number of audio files

    %% - - - - - - -  TFR computation - - - - - - - 
        % Computing the STFT:------------------------------------------------------
        tic_all = tic;
        [STFT,F,time] = spectrogram(x, hanning(Nf),Nf-hop,Nf,fs); 

        X_original = abs(STFT); 
        X = X_original;
        
        % Original FEMD v1 structure tensor:------------------------------------------

            % Using imgradientxy:
            [Xn, Xw] = imgradientxy(compress(X));

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
                    if (eig_val(1) + eig_val(4)) > epsilon
                        C(k,m) = ((eig_val(4) - eig_val(1)) / (eig_val(1) + eig_val(4)))^2;
                    end            
                end        
            end 
        
        % Using MRFCI v2 structure tensor:------------------------------------------
        %[C, v1, v2] = structure_tensor_v2(X,Ngauss,sigma_n,sigma_k);
        %alphas = atan(v2./v1);

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

            % Points: (how many points?)
            %[xpdf,ypdf] = fplot(pdfest, [-5,5], 'k', 100);
            fp = fplot(pdfest, [-10,10], 'k', 100);
            ypdf = fp.YData;
            xpdf = fp.XData;
            
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

        % Computing theoretical values of alpha:    
        if useTheoAlphas
            if NumOfSources == 1            
                alphasTEO = computeAlphasTEO(x1,fs,f1(1,:), (size(X,1)-1)*2, hop);
            elseif NumOfSources == 2
                alphasTEO = [computeAlphasTEO(x1, fs, f1(1,:), Nf, hop)  computeAlphasTEO(x2, fs, f2(1,:), Nf, hop)];
            end
            %alphasEST = alphasTEO;
        end

        % Compute STFChTs using alphas: - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        tic_STFT = tic;
        [X_2, F_2, time_2] = spectrogram(x, hanning(Nf_2), Nf_2 - hop_2, Nf_2, fs); 
        F_final = F_2;
        X_2 = abs(X_2);
        toc(tic_STFT)
        
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
            cat_spectrograms_interp = X_3;
            
            for ind = 1:size(cat_spectrograms,3)
                cat_spectrograms_interp = cat(3, cat_spectrograms_interp, interp2(m_grid_2, k_grid_2, ...
                                                            cat_spectrograms(:,:,ind), m_grid, k_grid, 'linear', 0));
            end
            F_final = F_3;
            cat_spectrograms = cat_spectrograms_interp;
        end
        
        cat_spectrograms = abs(cat_spectrograms( F_final <= ylim_vec(2)*1000, :, :)).^2;
    %     cat_spectrograms = abs(cat_spectrograms( F_final <= ylim_vec(2)*1000, :, :));
        F_final = F_final(F_final <= ylim_vec(2)*1000);
        
        if NumOfSources > 1
            switch comb_method
                case 1 % SWGM
                    STFChT_comb = SWGM_comb(cat_spectrograms, beta);                       
                case 2 % LS
                    STFChT_comb = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
                case 3 % SLS
                    STFChT_comb = spectrogram_comb_local_sparsity(cat_spectrograms, size_w_S, size_w_E, zeta, eta);
            end
        else
            STFChT_comb = STFChT.^2;
        end
    

end
