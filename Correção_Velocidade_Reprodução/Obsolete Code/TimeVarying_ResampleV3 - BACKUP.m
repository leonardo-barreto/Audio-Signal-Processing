%Funcao para reamostragem variante no tempo
%Entradas: vetor com fatores de reamostragem (resampleFactors), Sinal a ser reamostrado
%(inputSignal), número de coeficientes do filtro (filterCoeffs), tamanho da janela (frameSize)

%Saidas: Sinal reamostrado: (outputSignal)

%SIZE: MEDIDAS EM AMOSTRAS
%LENGTH: MEDIDAS EM TEMPO


function y = TimeVarying_ResampleV3(inputSignal,originalSamplingRate,resampleFactors,filterCoeffs,frameSize,DEBUG)
    
    %Debug and error checking
    
    DEBUG_OFF = 0;
    DEBUG_FULL = 1;

    %TOTAL BLOCKS?
    totalFrames = length(inputSignal)/frameSize;

    if (~CheckInteger(totalFrames))
        totalFrames = ceil(totalFrames);
        numberOfSamples = frameSize*totalFrames;
        inputSignal(1,length(inputSignal):numberOfSamples) = 0;
    end
    inputSignalSize = length(inputSignal); %Original signal size. Zero-padding at the end if needed to complete the final frame.

    if length(resampleFactors) ~= totalFrames
        X=sprintf('Error! Array of resampling factors must have the same size of the number of signal frames: %i.',totalFrames);
        error(X);
    end

    if DEBUG == DEBUG_FULL
        X = sprintf('Input Signal Size: %i', inputSignalSize);
        disp(X)
        X = sprintf('Frame Size: %i', frameSize);
        disp(X)
        X = sprintf('Total Frames: %i', totalFrames);
        disp(X)
        pause(1.0);
    end

    %Original signal parameters
    
    originalPeriod = 1/originalSamplingRate; %Sampling period of the original signal.
    frameIndex = 1; %Index of current frame (original signal).
    newSamplePosition = 1; %New sample position: moves proportional to time (normalized by the original period).
    resampledOutputIndex = 1; %Output signal index: moves without respect to time, just the overall samples.

    %Outputs

    outputSignal = zeros(1,2); %To prevent undefined variable errors.

    %Code

    while (frameIndex <= totalFrames)

        [stepSize,blockNewSamples,currentPeriod] = GenerateNormalizedParametersV2(resampleFactors,frameIndex,originalPeriod,blockSize,DEBUG);

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

        while (blockCurrentNewSample <= (blockNewSamples+1) & newSamplePosition <= inputSignalSize)
            
            if newSamplePosition > inputSignalSize 
                error('New sample position exceeds the limit.');
            end

            convolutionIndex = 1; %Internal index: moves from the center of sinc to borders to computate multiplications.
            sample = 0; %Current sample's value.

            if CheckInteger(newSamplePosition) % If newSamplePosition is an integer, it means it coincides with an original sample (time is normalized by originalPeriod!).
                outputSignal(resampledOutputIndex) = inputSignal(newSamplePosition);
            else
                leftClosestSample = floor(newSamplePosition);
                rightClosestSample = ceil(newSamplePosition);
                leftGap = newSamplePosition - leftClosestSample;
                rightGap = rightClosestSample - newSamplePosition;

                while convolutionIndex <= round(filterCoeffs/2)
                    convolutionIndex_left = leftClosestSample - (convolutionIndex-1); % From center to left
                    convolutionIndex_right = rightClosestSample + (convolutionIndex-1); % From center to right

                    if (convolutionIndex_left > 0 & convolutionIndex_left < inputSignalSize)
                        sample = sample + inputSignal(convolutionIndex_left)*sinc(-leftGap - (convolutionIndex-1));
                    end
                    if convolutionIndex_right <= inputSignalSize
                        sample = sample + inputSignal(convolutionIndex_right)*sinc(rightGap + (convolutionIndex-1));
                    end
                    convolutionIndex = convolutionIndex + 1;
                end
                outputSignal(resampledOutputIndex) = sample;
            end
            newSamplePosition = newSamplePosition + currentPeriod + (blockCurrentNewSample-1)*stepSize; %Interpolates between factors if needed.
            blockCurrentNewSample = blockCurrentNewSample + 1;
            resampledOutputIndex = resampledOutputIndex + 1;
        end
        frameIndex = frameIndex + 1;
    end
    y = outputSignal;
end