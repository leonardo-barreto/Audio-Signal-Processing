function [spectrg_SS, spectrg_Tr] = SSE_filter(spectrg, nFilter_SS, nFilter_Tr);

    % This function aims at producing Steady-State and Transient-enhanced spectrograms of a signal,
    % based on SSE filtering along both time and frequency dimensions.
    % 
    % Inputs:
    %   spectrg : input spectrogram.
    %   nFilter_SS : number of SSE filter coefficients for filtering along time axis (SS enhancement).
    %   nFilter_Tr : number of SSE filter coefficients for filtering along frequency axis (Transient enhancement).
    %
    % Outputs:
    %   spectr_SS : SS-enhanced version of input spectrogram.
    %   spectr_Tr : Transient-enhanced version of input spectrogram.
    %

    if mod(nFilter_SS,2) == 0 | mod(nFilter_Tr,2) == 0
        error('Filter lengths must be odd numbers.')
    end

    %% Time filtering (SS enhancement)

        % Extending the spectrogram by half a window size in order to mitigate border effects.
        mirrorLength = floor(nFilter_SS/2);
        startMirror = flip(spectrg(:,1:mirrorLength),2);
        %size(startMirror)
        endMirror = flip(spectrg(:,(end-mirrorLength)+1:end),2);
        %size(endMirror)
        spectrg_extended = [startMirror spectrg endMirror];
        %size(spectrg_extended)
        extendedLength = size(spectrg_extended,2); % length of extended spectrogram

        R_inverted = 1./(spectrg_extended);
        R_filtered = movmean(R_inverted,nFilter_SS,2);
        R_SS = 1./R_filtered;
        %size(R_SS)
        R_SS = R_SS(:,mirrorLength+1 : extendedLength-mirrorLength); %Discarding the mirrored edges
        %size(R_SS)
        
    %% Frequency filtering (Transient enhancement)

        % Extending the spectrogram by half a window size in order to mitigate border effects.
        mirrorLength = floor(nFilter_Tr/2);
        startMirror = flip(spectrg(1:mirrorLength,:),1);
        %size(startMirror)
        endMirror = flip(spectrg((end-mirrorLength)+1:end,:),1);
        %size(endMirror)
        spectrg_extended = [startMirror ; spectrg ; endMirror];
        %size(spectrg_extended)
        extendedLength = size(spectrg_extended,1); % length of extended spectrogram

        R_inverted = 1./spectrg_extended;
        R_filtered = movmean(R_inverted,nFilter_Tr,1);
        R_Tr = 1./R_filtered;
        %size(R_Tr)
        R_Tr = R_Tr(mirrorLength+1 : extendedLength-mirrorLength,:); %Discarding the mirrored edges
        %size(R_Tr)
    %% Wiener masking
        mask_SS=R_SS.^2./(R_SS.^2+R_Tr.^2);
        mask_Tr=R_Tr.^2./(R_SS.^2+R_Tr.^2);
        spectrg_SS=spectrg.*mask_SS;
        spectrg_Tr=spectrg.*mask_Tr;

end