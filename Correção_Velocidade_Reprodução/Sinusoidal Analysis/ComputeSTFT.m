function [spectrg,timeInstants] = ComputeSTFT(inputSignal,samplingRate,windowSize,overlapPerc,DEBUG)

    %overlapSize is given in SAMPLES.

    fftPoints = windowSize;

    overlapSize = (overlapPerc/100)*windowSize;

    [spectrg,normalizedFreqs,timeInstants] = spectrogram(inputSignal,hanning(windowSize),overlapSize,fftPoints,samplingRate);

    if DEBUG == 1;
        %spectrogram(inputSignal,windowSize,overlapSize,fftPoints,samplingRate,'yaxis');
        t = 0:max(timeInstants);
        f = 0:max(normalizedFreqs);
        surf(timeInstants,normalizedFreqs,abs(spectrg),'EdgeColor','none');   
        axis xy; 
        axis tight; 
        colormap('default'); 
        view(0,90);
    end

    %plot(timeInstants,)

end