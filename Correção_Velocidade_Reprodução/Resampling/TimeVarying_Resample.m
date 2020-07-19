%Funcao para reamostragem variante no tempo
%Entradas: vetor com fatores de reamostragem (resampleFactors), Sinal a ser reamostrado
%(inputSignal), número de coeficientes do filtro (filterCoeffs), parâmetros da análise senoidal (sinAnalysisParameters)

%Saidas: Sinal reamostrado: (outputSignal)

%SIZE: MEDIDAS EM AMOSTRAS
%LENGTH: MEDIDAS EM TEMPO


function outputSignal = TimeVarying_Resample(inputSignal,sinAnalysisParameters,resampleFactors,filterCoeffs,DEBUG)

    fprintf('\n\n------- TIME-VARYING RESAMPLING STARTED ------\n\n');
    
    %Debug and error checking
    
        DEBUG_OFF = 0;
        DEBUG_FULL = 1;
        PHASE_CORRECT = 0;

    %Gathering Sinusoidal Analysis data.

        originalSamplingRate = sinAnalysisParameters.samplingRate;
        timeInstants = sinAnalysisParameters.timeInstants;
        frameSize = sinAnalysisParameters.windowSize;
        totalFrames = sinAnalysisParameters.totalFrames;
        hopSize = sinAnalysisParameters.hopSize;

    %Original signal parameters
    
        inputSignalSize = length(inputSignal);
        originalPeriod = 1/originalSamplingRate; %Sampling period of the original signal.

    if DEBUG == DEBUG_FULL
        X = sprintf('Input Signal Size: %i', inputSignalSize);
        disp(X)
        X = sprintf('Frame Size: %i', frameSize);
        disp(X)
        X = sprintf('Total Frames: %i', totalFrames);
        disp(X)
        pause(1.0);
    end

    %Building the center samples of each block: the first and last frames have to be accounted for.
        centerSamples = [1,timeInstants*originalSamplingRate,inputSignalSize];

    %Outputs
        outputSignal = []; %To prevent undefined variable errors.

    %Phase Correction
        if PHASE_CORRECT == 1
            phaseCorrection = 0; %Factor that will correct the synchrony between processed blocks.
        end

    %Resampling

        frameIndex = 1; %Index of current frame (original signal).
        newSamplePosition = 1; %New sample position: moves proportional to time (normalized by the original period).
        resampledOutputIndex = 1; %Output signal index: moves without respect to time, just the overall samples.

        while (frameIndex <= totalFrames)

            blockSize = centerSamples(frameIndex+1)-centerSamples(frameIndex)-1;

            [stepSize,blockNewSamples,currentPeriod] = GenerateNormalizedParameters(resampleFactors,frameIndex,originalPeriod,blockSize,DEBUG);

            if DEBUG == DEBUG_FULL
                X = sprintf('Current Frame: %i of %i', frameIndex, totalFrames);
                disp(X)
                X = sprintf('Step Size: %i', stepSize);
                disp(X)
                X = sprintf('New Samples: %i', blockNewSamples);
                disp(X)
                pause(0.5);
            end 

            blockCurrentNewSample = 1;

            if newSamplePosition > inputSignalSize 
                error('New sample position exceeds the limit (greater than signal size).');
            end

            if (newSamplePosition < 1)
                error('New sample position exceeds the limit (less than signal starting point)');
            end

            while (blockCurrentNewSample <= (blockNewSamples+1) & newSamplePosition <= inputSignalSize)
                
                convolutionIndex = 0; %Internal index: moves from the center of sinc to borders to compute multiplications.
                sample = 0; %Current sample's value.

                if CheckInteger(newSamplePosition) % If newSamplePosition is an integer, it means it coincides with an original sample (time is normalized by originalPeriod!).
                    outputSignal(resampledOutputIndex) = inputSignal(newSamplePosition);
                else
                    leftClosestSample = floor(newSamplePosition);
                    rightClosestSample = ceil(newSamplePosition);
                    leftGap = newSamplePosition - leftClosestSample;
                    rightGap = rightClosestSample - newSamplePosition;

                    while convolutionIndex <= round(filterCoeffs/2)
                        convolutionIndex_left = leftClosestSample - convolutionIndex; % From center to left
                        convolutionIndex_right = rightClosestSample + convolutionIndex; % From center to right

                        if (convolutionIndex_left >= 1 & convolutionIndex_left < inputSignalSize)
                            sample = sample + inputSignal(convolutionIndex_left)*sinc(-leftGap - convolutionIndex); %MATLAB's sinc zeroes itself each step of 1, so we can work in samples directly.
                        end
                        if convolutionIndex_right <= inputSignalSize
                            sample = sample + inputSignal(convolutionIndex_right)*sinc(rightGap + convolutionIndex);
                        end
                        convolutionIndex = convolutionIndex + 1;
                    end
                    outputSignal(resampledOutputIndex) = sample;
                end
                newSamplePosition = newSamplePosition + currentPeriod + (blockCurrentNewSample-1)*stepSize; %Interpolates between factors if needed.
                blockCurrentNewSample = blockCurrentNewSample + 1;
                resampledOutputIndex = resampledOutputIndex + 1;
            end
            if (PHASE_CORRECT == 1)
                phaseCorrection = blockSize+1 - (blockNewSamples+1)*currentPeriod - (blockNewSamples+1)*stepSize*(blockNewSamples/2)
                phaseCorrectionMY = (blockSize+1) - blockNewSamples*(currentPeriod + stepSize*(blockNewSamples-1)/2)/originalPeriod
                newSamplePosition = newSamplePosition - phaseCorrection
            end
            frameIndex = frameIndex + 1;
        end

        fprintf('\n\n------- TIME-VARYING RESAMPLING FINISHED ------\n\n');

end