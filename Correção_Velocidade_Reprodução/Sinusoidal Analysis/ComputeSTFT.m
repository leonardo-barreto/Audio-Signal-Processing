function [spectrgMatrix,freqComponents,frameTimeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints)

    % This function computes a signal's STFT (by FFT).

    DEBUG = 0;
    METHOD = 0;
    
    if METHOD == 0
        disp(' Method used: MATLAB spectrogram.');
    else
        disp(' Method used: stft function from MathWorks site.');
    end

    switch windowType
        case 'hann' 
            windowFunction = hann(windowSize,'periodic');
        case 'hamming'
            windowFunction = hamming(windowSize,'periodic');
        otherwise
            error('Invalid window type. Valid options are ''hann'' or ''hamming''.');
    end

    if METHOD == 0
        overlapSize = floor((overlapPerc/100)*windowSize); %This converts the overlap percentage to actual overlap size.
        [s,f,t,ps] = spectrogram (inputSignal,windowFunction,overlapSize,2*fftPoints-1,'power','onesided');
        freqComponents = (samplingRate/2*pi).*f;
        powerMatrix = ps;
    else
        hopSize = floor(((100-overlapPerc)/100)*windowSize); %This converts the overlap percentage to hop size.
        [s, f, t] = stft(inputSignal,windowFunction,hopSize,fftPoints,samplingRate);
        freqComponents = f;
        powerMatrix = power(abs(s),2)/fftPoints;
    end

    spectrgMatrix = s;
    frameTimeInstants = t;

    if DEBUG == 1
        % plot the spectrogram
        S = 20*log10(powerMatrix);
        figure(1)
        surf(t, f, S)
        shading interp
        axis tight
        view(0, 90)
        set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
        xlabel('Time, s')
        ylabel('Frequency, Hz')
        title('Amplitude spectrogram of the signal')

        hcol = colorbar;
        set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
        ylabel(hcol, 'Magnitude, dB')
    end
    

end