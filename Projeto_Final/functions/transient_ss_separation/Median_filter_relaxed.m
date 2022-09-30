function [ spectrg_SS, spectrg_Tr ] = Median_filter_relaxed(spectrg,nFilter_SS,nFilter_Tr)

    % This function aims at producing Steady-State and Transient-enhanced spectrograms of a signal,
    % based on median filtering along both time and frequency dimensions.
    %  
    % Original code implementaion made by Ignacio Irigaray 
    % (Universidad de la Republica, Montevideo, Uruguay)
    % 
    % Inputs:
    %   spectrg    : Spectrogram module matrix
    %   nFilter_SS : filter length of steady-state median filter 
    %   nFilter_Tr : filter length of transient median filter
    % 
    % Outputs:
    %   spectrg_SS : Steady-state components spectrogram
    %   spectrg_Im : Transient components spectrogram
    %

    if mod(nFilter_SS,2) == 0 | mod(nFilter_Tr,2) == 0
        nFilter_SS = nFilter_SS + 1;
        nFilter_Tr = nFilter_Tr + 1;
        warning('Filter lengths must be odd. Using %i and %i instead.',nFilter_SS,nFilter_Tr)
    end

    nFilter_SS = (nFilter_SS-1)/2; % Total length to neighbourhood size
    nFilter_Tr = (nFilter_Tr-1)/2;

    spectrg_SS=zeros(size(spectrg));
    spectrg_Tr=zeros(size(spectrg));

    timeLength = size(spectrg,2);
    freqHeight = size(spectrg,1);

    %% Median filtering along time
        ext = floor(nFilter_SS/4)+1;

        for  k = 1:timeLength
            if (k > nFilter_SS) && (k <= timeLength-nFilter_SS)
                for i = 1:freqHeight
                    if (i > ext) && (i <= freqHeight-ext)
                        spectrg_SS(i,k) = median(max(spectrg(i-ext:i+ext,k-nFilter_SS:k+nFilter_SS)),2);
                    else
                        spectrg_SS(i,k) = median(spectrg(i,k-nFilter_SS:k+nFilter_SS),2);
                    end
                end
            else
                if k <= nFilter_SS
                    for i = 1:freqHeight
                        if (i > ext) && (i <= freqHeight-ext)
                            spectrg_SS(i,k) = median(max(spectrg(i-ext:i+ext,1:k+nFilter_SS)),2);
                        else
                            spectrg_SS(i,k) = median(spectrg(i,1:k+nFilter_SS),2);
                        end
                    end
                else
                    for i = 1:freqHeight
                        if (i > ext) && (i <= freqHeight-ext)
                            spectrg_SS(i,k) = median(max(spectrg(i-ext:i+ext,k-nFilter_SS:k)),2);
                        else
                            spectrg_SS(i,k) = median(spectrg(i,k-nFilter_SS:k),2);
                        end
                    end
                end
                
            end
        end

    %% Median filtering along frequency
        disp('done with time')
        for k = 1:freqHeight
            if (k > nFilter_Tr) && (k <= freqHeight-nFilter_Tr)
                spectrg_Tr(k,:) = median(spectrg(k-nFilter_Tr:k+nFilter_Tr,:),1);
            else
                if k <= nFilter_Tr
                    spectrg_Tr(k,:) = median(spectrg(1:k+nFilter_Tr,:),1);
                else    
                    spectrg_Tr(k,:) = median(spectrg(k-nFilter_Tr:k,:),1);
                end
                %spectrg_Tr(k,:)=0;
            end
        end


    %% Wiener masks
        mask_SS=spectrg_SS.^2./(spectrg_SS.^2+spectrg_Tr.^2);
        mask_Tr=spectrg_Tr.^2./(spectrg_SS.^2+spectrg_Tr.^2);
        spectrg_SS=spectrg.*mask_SS;
        spectrg_Tr=spectrg.*mask_Tr;

end

