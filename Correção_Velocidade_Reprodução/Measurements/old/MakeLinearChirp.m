function outputSignal = MakeLinearChirp(initialFreq,finalFreq,sweepStart,sweepEnd,totalTime,samplingRate)

    sweepTime = sweepEnd-sweepStart;
    
    timeAxis0 = 0:(1/samplingRate):sweepStart-1/samplingRate;
    timeAxisChirp = sweepStart:(1/samplingRate):sweepEnd-1/samplingRate;
    timeAxis1 = sweepEnd:1/samplingRate:totalTime-1/samplingRate;

    chirpyness = (finalFreq-initialFreq)/sweepTime;

    waveform0 = sin(2*pi*initialFreq*timeAxis0);
    waveformChirp = sin(2*pi*((chirpyness/2)*(timeAxisChirp.^2)+initialFreq*timeAxisChirp));
    waveform1 = sin(2*pi*finalFreq*timeAxis1);

    outputSignal = [waveform0 waveformChirp waveform1]; 