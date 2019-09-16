%Funcao para reamostragem variante no tempo
%Entradas: vetor com fatores de reamostragem (resampleFactors), Sinal a ser reamostrado
%(inputSignal), número de coeficientes do filtro (filterCoeffs), tamanho da janela (windowSize)

%Saidas: Sinal reamostrado: (outputSignal)

%SIZE: MEDIDAS EM AMOSTRAS
%LENGTH: MEDIDAS EM TEMPO


function y = TimeVarying_ResampleV3(inputSignal,originalSamplingRate,resampleFactors,filterCoeffs,windowSize,DEBUG)
    
    %Debug and error checking
    
    DEBUG_OFF = 0;
    DEBUG_FULL = 1;

    totalWindows = length(inputSignal)/windowSize;
    if (~CheckInteger(totalWindows))
        totalWindows = ceil(totalWindows);
        numberOfSamples = windowSize*totalWindows;
        inputSignal(1,length(inputSignal):numberOfSamples) = 0;
    end
    inputSignalSize = length(inputSignal); %Original signal size. Zero-padding at the end if needed to complete the final window.

    if length(resampleFactors) ~= totalWindows
        X=sprintf('Error! Array of resampling factors must have the same size of the number of signal blocks: %i.',totalWindows);
        error(X);
    end

    if DEBUG == DEBUG_FULL
        X = sprintf('Input Signal Size: %i', inputSignalSize);
        disp(X)
        X = sprintf('Window Size: %i', windowSize);
        disp(X)
        X = sprintf('Total Windows: %i', totalWindows);
        disp(X)
        pause(1.0);
    end

    %Original signal parameters
    
    blockSize = windowSize; %TEMP
    originalPeriod = 1/originalSamplingRate; %Sampling period of the original signal.
    windowIndex = 1; %Index of current window (original signal).
    newSamplePosition = 1; %New sample position: moves proportional to time (normalized by the original period).
    resampledOutputIndex = 1; %Output signal index: moves without respect to time, just the overall samples.

    %Outputs

    outputSignal = zeros(1,2); %To prevent undefined variable errors.

    %Code

    while (windowIndex <= totalWindows)

        [stepSize,blockNewSamples,currentPeriod] = GenerateNormalizedParametersV2(resampleFactors,windowIndex,originalPeriod,blockSize,DEBUG);

        if DEBUG == DEBUG_FULL
            X = sprintf('Current Window: %i of %i', windowIndex, totalWindows);
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
        windowIndex = windowIndex + 1;
    end
    y = outputSignal;
end