function [spectrgMatrix,powerSpecMatrix,freqComponents,timeInstants] = ComputeSTFT(inputSignal,samplingRate,windowSize,overlapPerc,DEBUG)

    fftPoints = windowSize; %DFT is calculated through all points of the window.

    overlapSize = (overlapPerc/100)*windowSize; %This converts the overlap percentage to actual overlap size.

    [spectrgMatrix,powerSpecMatrix,freqComponents,timeInstants] = spectrogram (inputSignal,hann(windowSize),overlapSize,fftPoints,samplingRate,'power','centered');

    if DEBUG == 1;
        spectrogram(inputSignal,windowSize,overlapSize,fftPoints,samplingRate,'yaxis');
        %t = 0:max(timeInstants);
        %f = 0:max(normalizedFreqs);
        %surf(timeInstants,normalizedFreqs,abs(spectrgMatrix),'EdgeColor','none');   
        %axis xy; 
        %axis tight; 
        %colormap('default'); 
        %view(0,90);
    end

end