function [ spectrg_SS, spectrg_Tr, spectrg_Res ] = Iterative_HPR_Separation(spectrg, nFilter_SS, nFilter_Tr, niter, method, varargin)

    % This function aims at producing Steady-State, Transient and Residual-enhanced spectrograms of a signal,
    % based on median filtering or SSE filtering along both time and frequency dimensions.
    %  
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
    Median_functions = {@Median_filter @Median_filter_relaxed};
    SSE_functions = {@SSE_filter @SSE_filter_relaxed};
    func_idx = 1;

    spectrg_SS = zeros(size(spectrg));
    spectrg_Tr = zeros(size(spectrg));
    spectrg_Res = zeros(size(spectrg));
    spectrg_SSTemp = zeros(size(spectrg));
    spectrg_TrTemp = zeros(size(spectrg));

    if (niter < 1)
        error('niter (number of iterations) must be greater than or equal to 1.');
    end
    if (nargin > 5)
        if (nargin > 6 && (~nnz(ismember(varargin,'dB')) | ~nnz(ismember(varargin,'relaxed'))))
            error('Invalid optional arguments inserted. Options are ''dB'' and/or ''relaxed''.');
        elseif (nnz(ismember(varargin,'dB')) && ~strcmpi(method,'SSE')) 
            error('Cannot have dB with non-SSE method.');
        elseif nnz(ismember(varargin,'relaxed')) %Relaxed components
            func_idx = 2;
        end
        if (nargin > 7)
        error('Max number of arguments including optionals is 7.');
        end
    end

    if (strcmpi(method,'median')) % Median method

        Median_function = Median_functions{func_idx};
        [spectrg_SS,spectrg_Tr] = Median_function(spectrg,nFilter_SS,nFilter_Tr); % First processing

        if (niter > 1) % Iterative processing
            for i = 2:niter
                [spectrg_SS,spectrg_TrTemp] = Median_function(spectrg_SS,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_TrTemp;

                [spectrg_SSTemp,spectrg_Tr] = Median_function(spectrg_Tr,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_SSTemp;
            end
        end

    elseif (strcmpi(method,'SSE')) % SSE method

        SSE_function = SSE_functions{func_idx};

        if nnz(ismember(varargin,'dB'))
            spectrg = spectrg/10;
            spectrg = power(10,spectrg);
        end

        [spectrg_SS,spectrg_Tr] = SSE_function(spectrg,nFilter_SS,nFilter_Tr); % First processing
        
        if (niter > 1) % Iterative processing
            for i = 2:niter
                [spectrg_SS,spectrg_TrTemp] = SSE_function(spectrg_SS,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_TrTemp;

                [spectrg_SSTemp,spectrg_Tr] = SSE_function(spectrg_Tr,nFilter_SS,nFilter_Tr);
                spectrg_Res = spectrg_Res + spectrg_SSTemp;
            end
        end
    else
        error('Method incorrect. It must be either ''median'' or ''SSE''.');
    end
    

end


