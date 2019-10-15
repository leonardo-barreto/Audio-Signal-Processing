function plot_fft(signal,Fs)

    Ts = 1/Fs;         
    N = length(signal);                                     
    signalLength = (N-1)*Ts;

    S = fft(signal);

    amplitude_spectrum = abs(S);

    frequency_domain = -Fs/2:(Fs/N):Fs/N;

    length(frequency_domain)
    length(amplitude_spectrum)

    if length(frequency_domain) ~= length(amplitude_spectrum)
        error('Error! Frequency domain must be the same size of the spectrum!');
    end

    figure
    plot(frequency_domain,amplitude_spectrum)
    title('Amplitude Spectrum');
    xlabel('Hz');
    ylabel('Amplitude');
end




    
 