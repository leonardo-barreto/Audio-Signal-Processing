%Funcao para reamostragem variante no tempo
%Entradas: vetor com fatores de reamostragem (resampleFactors), Sinal a ser reamostrado
%(inputSignal), número de coeficientes do filtro (filterCoeffs), tamanho do bloco de influência (blockSize)

%Saidas: Sinal reamostrado: (outputSignal)

%SIZE: MEDIDAS EM AMOSTRAS
%LENGTH: MEDIDAS EM TEMPO


function y = TimeVarying_ResampleV2(inputSignal,originalSamplingRate,resampleFactors,filterCoeffs,windowSize,DEBUG)
    
    %Definitions
    
    DEBUG_OFF = 0;
    DEBUG_ERRORS = 1;
    DEBUG_FULL = 2;

    windowCenter = (windowSize/2) + 1;
    
    blockSize = windowSize;

    totalWindows = length(inputSignal)/windowSize;
    if (~CheckInteger(totalWindows))
        totalWindows = ceil(totalWindows);
        numberOfSamples = windowSize*totalWindows;
        inputSignal(1,length(inputSignal):numberOfSamples) = 0;
    end
    inputSignalSize = length(inputSignal); %Tamanho do sinal de entrada com zero-padding ao final se necessário.

    if DEBUG == DEBUG_FULL
        X = sprintf('Input Signal Size: %i', inputSignalSize);
        disp(X)
        X = sprintf('Window Size: %i', windowSize);
        disp(X)
        X = sprintf('Total Windows: %i', totalWindows);
        disp(X)
        pause(0.5);
    end
    
    inputPeriod = 1/originalSamplingRate; %período de amostragem original
    windowIndex = 1; %índice que se move a cada janela processada.

    newSamplePosition = 1; %posição da nova amostra a ser posicionada: se move no tempo.

    resampledOutputIndex = 1; %índice do sinal de saída, que se move em amostras.

    %Outputs
    outputSignal = zeros(1,2); % Para prevenir erros de undefined variable

    if length(resampleFactors) ~= totalWindows
        X=sprintf('Error! Array of resampling factors must have the same size of the number of signal blocks: %i.',totalWindows);
        error(X);
    end


    while (windowIndex <= totalWindows)

        [stepSize,blockNewSamples,currentPeriod] = GenerateNormalizedParameters(resampleFactors,windowIndex,inputPeriod,blockSize);

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

        while (blockCurrentNewSample <= blockNewSamples & newSamplePosition <= inputSignalSize) % percorre um bloco de influência, que terá blockNewSamples novas amostras.
            
            if newSamplePosition > inputSignalSize 
                error('New sample position exceeds the limit.');
            end

            convolutionIndex = 1; % índice interno que percorre o sinal de entrada para a multiplicação pela sinc
            sample = 0;

            if CheckInteger(newSamplePosition) % se a amostra na nova taxa coincide com a amostra na taxa antiga
                outputSignal(resampledOutputIndex) = inputSignal(newSamplePosition);
            else
                leftClosestSample = floor(newSamplePosition);
                rightClosestSample = ceil(newSamplePosition);
                leftGap = newSamplePosition - leftClosestSample;
                rightGap = rightClosestSample - newSamplePosition;

                while convolutionIndex <= round(filterCoeffs/2)
                    convolutionIndex_left = leftClosestSample - (convolutionIndex-1); % do centro da sinc para esq
                    convolutionIndex_right = rightClosestSample + (convolutionIndex-1); % do centro da sinc para dir

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
            newSamplePosition = newSamplePosition + currentPeriod + (blockCurrentNewSample-1)*stepSize;
            blockCurrentNewSample = blockCurrentNewSample + 1;
            resampledOutputIndex = resampledOutputIndex + 1;
        end
        windowIndex = windowIndex + 1;
    end

    y = outputSignal;
end