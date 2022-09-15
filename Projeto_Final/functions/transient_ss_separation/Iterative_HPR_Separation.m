function [ spectrg_SS, spectrg_Tr, spectrg_Res ] = Iterative_HPR_Separation(varargin)%spectrg, nFilter_SS, nFilter_Tr, niter, method, varargin)

    % This function aims at producing Steady-State, Transient and Residual-enhanced spectrograms of a signal,
    % based on median filtering or SSE filtering along both time and frequency dimensions.
    %  
    % Base code implementation made by Ignacio Irigaray 
    % (Universidad de la Republica, Montevideo, Uruguay)
    %
    % Inputs:
    %   spectrg    : Spectrogram module matrix
    %   nFilter_SS : filter length of steady-state filters 
    %   nFilter_Tr : filter length of transient filters
    %   niter      : number of separation iterations
    %   method     : filtering method ('median' or 'SSE')
    %
    % Optional input arguments (varargin):
    %   'dB' (only necessary for SSE)  : indicates an input spectrogram in dB, which prompts a conversion back to non-dB
    %   'relaxed'   : indicates the use of Relaxed Components (allows some level of frequency variation in steady-state)
    % 
    % Outputs:
    %   spectrg_SS  : Steady-state components spectrogram
    %   spectrg_Im  : Transient components spectrogram
    %   spectrg_Res : Residual components spectrogram
    %

    %%Initializing variables
    spectrg_SS = zeros(size(spectrg));
    spectrg_Tr = zeros(size(spectrg));
    spectrg_Res = zeros(size(spectrg));
    spectrg_SSTemp = zeros(size(spectrg));
    spectrg_TrTemp = zeros(size(spectrg));

    if nnz(ismember(varargin,'relaxed'))
        %Relaxed components
    else
        error('You are writing the optional arguments wrong. Options are ''dB'' and ''relaxed''.');
    end

    if (strcmpi(method,'median')) % Median method

        [spectrg_SS,spectrg_Tr] = Median_filter(spectrg,nFilter_SS,nFilter_Tr); % First processing

        if (niter ~= 0) % Iterative processing
            for i = 1:niter
                [spectrg_SS,spectrg_TrTemp] = Median_filter(spectrg_SS,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_TrTemp;

                [spectrg_SSTemp,spectrg_Tr] = Median_filter(spectrg_Tr,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_SSTemp;
            end
        end

    elseif (strcmpi(method,'SSE')) % SSE method

        if nargin > 5 && nnz(ismember(varargin,'dB'))
            spectrg = spectrg/10;
            spectrg = power(10,spectrg);
        else 
            if nargin > 6
                error('You are writing the optional arguments wrong. Options are ''dB'' and ''relaxed''.');
            end
        end

        [spectrg_SS,spectrg_Tr] = SSE_filter(spectrg,nFilter_SS,nFilter_Tr); % First processing
        
        if (niter ~= 0) % Iterative processing
            for i = 1:niter
                [spectrg_SS,spectrg_TrTemp] = SSE_filter(spectrg_SS,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_TrTemp;

                [spectrg_SSTemp,spectrg_Tr] = SSE_filter(spectrg_Tr,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_SSTemp;
            end
        end
    else
        error('Method incorrect. It must be either ''median'' or ''SSE''.');
    end
    

end


