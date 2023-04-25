function [TFR_CQT, f, t] = TFAnalysis_CQT(x,fs)

    %   This function makes a time-frequency analysis of a signal, using the Constant-Q Transform.

    %% - - - - - - Parameters configuration - - - - - - 

        % Analysis TFR - - - - - - - - - - - - - - - - - - - -

        hop = 128;
        N_w = 4096;
        NFFT = N_w;
        bins_per_octave = 24;
        fmax = fs/2;     %center frequency of the highest frequency bin 
        fmin = fmax/512; %lower boundary for CQT (lowest frequency bin will be immediately above this): fmax/<power of two> 


    %% - - - - - - -  TFR computation - - - - - - -
        Xcqt = cqtPerfectRast(x,fmin,fmax,bins_per_octave,fs,'q',1,'atomHopFactor',0.25,'thresh',0.0005,'win','sqrt_blackmanharris');

        emptyHops = Xcqt.intParams.firstcenter/Xcqt.intParams.atomHOP;
        maxDrop = emptyHops*2^(Xcqt.octaveNr-1)-emptyHops;
        droppedSamples = (maxDrop-1)*Xcqt.intParams.atomHOP + Xcqt.intParams.firstcenter;
        
        absCQT = abs(Xcqt.spCQT);
        t = (1:size(absCQT,2))*Xcqt.intParams.atomHOP-Xcqt.intParams.preZeros+droppedSamples;
        t = t./fs;
        f = 1:size(absCQT,1); 
        spectrgMatrix = absCQT;
        
        TFR_CQT = power(abs(spectrgMatrix),2);%/NFFT;
        %TFR_CQT = abs(spectrgMatrix);%/NFFT;
    
end