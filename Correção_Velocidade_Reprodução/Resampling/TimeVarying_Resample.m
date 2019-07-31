%Funcao para reamostragem variante no tempo
%Entradas: vetor com fatores de reamostragem (resampleFactors), Sinal a ser reamostrado
%(inputSignal), número de coeficientes do filtro (filterCoeffs), tamanho do bloco de influência (blockSize)

%Saidas: Sinal reamostrado: (outputSignal)

%SIZE: MEDIDAS EM AMOSTRAS
%LENGTH: MEDIDAS EM TEMPO


function y = TimeVarying_Resample(inputSignal,originalSamplingRate,resampleFactors,filterCoeffs,blockSize,DEBUG)
    
    %Definitions
    
    DEBUG_OFF = 0;
    DEBUG_ERRORS = 1;
    DEBUG_FULL = 2;

    numberOfBlocks = length(inputSignal)/blockSize;
    if (~CheckInteger(numberOfBlocks))
        numberOfBlocks = ceil(numberOfBlocks);
        numberOfSamples = blockSize*numberOfBlocks;
        inputSignal(1,length(inputSignal):numberOfSamples) = 0;
    end
    inputSignalSize = length(inputSignal); %Tamanho do sinal de entrada com zero-padding ao final se necessário.

    if DEBUG == DEBUG_FULL
        X = sprintf('Input Signal Size: %i', inputSignalSize);
        disp(X)
        pause(0.5);
    end
    
    inputPeriod = 1/originalSamplingRate; %período de amostragem original
    windowIndex = 1; %índice que se move a cada janela processada.

    newSamplePosition = 1; %posição da nova amostra a ser posicionada: se move no tempo.

    resampledIndex = 1; %índice do sinal de saída, que se move em todas as diferentes taxas.

    %Outputs
    outputSignal = zeros(1,2); % Para prevenir erros de undefined variable

    if length(resampleFactors) ~= numberOfBlocks
        X=sprintf('Error! Array of resampling factors must have the same size of the number of signal blocks: %i.',numberOfBlocks);
        error(X);
    end


    while windowIndex <= numberOfBlocks

        [stepSize,numberOfSteps,currentPeriod] = GenerateNormalizedParameters(resampleFactors,windowIndex,inputPeriod,blockSize);

        currentStep = 0;

        while (currentStep < numberOfSteps & newSamplePosition <= inputSignalSize) % percorre um bloco de influência
            
            if newSamplePosition > inputSignalSize 
                error('New sample position exceeds the limit.');
            end

            convolutionIndex = 1; % índice interno que percorre o sinal de entrada para a multiplicação pela sinc
            sample = 0;
            if CheckInteger(newSamplePosition) % se a amostra na nova taxa coincide com a amostra na taxa antiga
                outputSignal(resampledIndex) = inputSignal(newSamplePosition);
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
                outputSignal(resampledIndex) = sample;
            end
            newSamplePosition = newSamplePosition + currentPeriod + currentStep*stepSize;
            currentStep = currentStep + 1;
            resampledIndex = resampledIndex + 1;
        end
        windowIndex = windowIndex + 1;
    end

    y = outputSignal;
end