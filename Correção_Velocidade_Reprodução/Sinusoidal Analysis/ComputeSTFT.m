function [spectrgMatrix,freqComponents,frameTimeInstants,powerMatrix] = ComputeSTFT(inputSignal,samplingRate,windowType,windowSize,overlapPerc,fftPoints)

    % This function computes a signal's STFT (by FFT).

    DEBUG = 1;
    METHOD = 1;

    switch windowType
        case 'hann' 
            windowFunction = hann(windowSize,'periodic');
        case 'hamming'
            windowFunction = hamming(windowSize,'periodic');
        otherwise
            error('Invalid window type. Valid options are ''hann'' or ''hamming''.');
    end

    if METHOD == 0
        fprintf(' Method used: MATLAB spectrogram.\n');
        overlapSize = floor((overlapPerc/100)*windowSize); %This converts the overlap percentage to actual overlap size.
        [s,f,t,ps] = spectrogram (inputSignal,windowFunction,overlapSize,fftPoints,'power','onesided');
        freqComponents = transpose((samplingRate/2*pi).*f);
        powerMatrix = ps;
    else
        fprintf(' Method used: stft function from Hristo Zhivomirov (MathWorks site).\n');
        hopSize = floor(((100-overlapPerc)/100)*windowSize); %This converts the overlap percentage to hop size.
        [s, f, t] = stft(inputSignal,windowFunction,hopSize,fftPoints,samplingRate);
        freqComponents = f;
        powerMatrix = power(abs(s),2)/fftPoints;
    end

    spectrgMatrix = s;
    frameTimeInstants = t;

    if DEBUG == 1
        % plot the spectrogram
        S = 10*log10(powerMatrix);
        figure(1)
        surf(t, f, S)
        shading interp
        axis tight
        view(0, 90)
        set(gca, 'FontSize', 30,'yscale','log')
        xlabel('Tempo (s)','FontSize', 30)
        ylabel('Frequencia (Hz)','FontSize', 30)
        title('Espectrograma de potencia')

        hcol = colorbar;
        set(hcol, 'FontSize', 30)
        ylabel(hcol, 'Potencia (dB)')
        %ylim([0 3000])
    end
    

end