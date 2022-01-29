function [spectrg_SS, spectrg_Tr] = SSE_filter(spectrg, nSSE_SS, nSSE_Tr, dBTrue);

    % This function aims at producing Steady-State and Transient-enhanced spectrograms of a signal,
    % based on SSE filtering along both time and frequency dimensions.
    % 
    % Inputs:
    %   spectrg : input spectrogram.
    %   nSSE_SS : number of SSE filter coefficients for filtering along time axis (SS enhancement).
    %   nSSE_Tr : number of SSE filter coefficients for filtering along frequency axis (Transient enhancement).
    %   dBTrue : 1 if input spectrogram is in dB, otherwise 0.
    %
    % Outputs:
    %   spectr_SS : SS-enhanced version of input spectrogram.
    %   spectr_Tr : Transient-enhanced version of input spectrogram.
    %

    % Gathering frame data
    if dBTrue == 1
        spectrg = spectrg/10;
        spectrg = power(10,spectrg);
    end

    %% Time filtering (SS enhancement)

        % Extending the spectrogram by half a window size in order to mitigate border effects.
        mirrorLength = ceil(nSSE_SS/2);
        startMirror = flip(spectrg(:,1:mirrorLength),2);
        endMirror = flip(spectrg(:,(end-mirrorLength)+1:end),2);
        spectrg_extended = [startMirror spectrg endMirror];
        extendedLength = size(spectrg_extended,2); % length of extended spectrogram

        R_inverted = 1./spectrg_extended;
        R_filtered = movmean(R_inverted,nSSE_SS,2);
        R_SS = 1./R_filtered;
        R_SS = R_SS(:,length(startMirror)+1 : extendedLength-length(endMirror)); %Discarding the mirrored edges

    %% Frequency filtering (Transient enhancement)

        % Extending the spectrogram by half a window size in order to mitigate border effects.
        mirrorLength = ceil(nSSE_Tr/2);
        startMirror = flip(spectrg(1:mirrorLength,:),1);
        endMirror = flip(spectrg((end-mirrorLength)+1:end,:),1);
        spectrg_extended = [startMirror ; spectrg ; endMirror];
        extendedLength = size(spectrg_extended,1); % length of extended spectrogram

        R_inverted = 1./spectrg_extended;
        R_filtered = movmean(R_inverted,nSSE_Tr,1);
        R_Tr = 1./R_filtered;
        R_Tr = R_Tr(length(startMirror)+1 : extendedLength-length(endMirror),:); %Discarding the mirrored edges

    %% Wiener masking
        mask_SS=R_SS.^2./(R_SS.^2+R_Tr.^2);
        mask_Tr=R_Tr.^2./(R_SS.^2+R_Tr.^2);
        spectrg_SS=spectrg.*mask_SS;
        spectrg_Tr=spectrg.*mask_Tr;
    

end