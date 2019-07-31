%Funcao para reamostragem variante no tempo
%Entradas: vetor com fatores de reamostragem (resampleFactors), Sinal a ser reamostrado
%(inputSignal), número de coeficientes do filtro (filterCoeffs), tamanho do bloco de influência (blockSize)

%Saidas: Sinal reamostrado: (outputSignal)

%SIZE: MEDIDAS EM AMOSTRAS
%LENGTH: MEDIDAS EM TEMPO


function y = ResampleSignal_TimeVarying(inputSignal,originalSamplingRate,resampleFactors,filterCoeffs,blockSize)
    
    %Definitions
    inputSignalSize = length(inputSignal);
    inputPeriod = 1/originalSamplingRate; %período de amostragem original
    windowIndex = 1; %índice que se move a cada janela processada.

    newSamplePosition = 1; %posição da nova amostra a ser posicionada: se move no tempo.

    resampledIndex = 1; %índice do sinal de saída, que se move em todas as diferentes taxas.
    
    blockLength = blockSize*inputPeriod; %tamanho do bloco de influência de cada fator.

    %Outputs
    outputSignal = zeros(1,1); % Para prevenir erros de undefined variable


    while windowIndex <= length(resampleFactors)
        
        currentFactor = resampleFactors(windowIndex);
        nextFactor = resampleFactors(windowIndex+1);

        stepSize = ((nextFactor^2) - (currentFactor^2))/(2*blockLength-(currentFactor+nextFactor));
        if stepSize > 0
            numberOfSteps = floor(2*blockLength/(currentFactor+nextFactor));
        else
            numberOfSteps = ceil(2*blockLength/(currentFactor+nextFactor));
        end

        currentPeriod = 1/currentFactor;
        currentStep = 1;

        while currentStep <= numberOfSteps % percorre um bloco de influência
            newSamplePosition = newSamplePosition + currentPeriod + (currentStep-1)*stepSize;
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
                    if convolutionIndex_left > 0
                        sample = sample + inputSignal(convolutionIndex_left)*sinc(-leftGap - (convolutionIndex-1));
                    end
                    if convolutionIndex_right < inputSignalSize
                        sample = sample + inputSignal(convolutionIndex_right)*sinc(rightGap + (convolutionIndex-1));
                    end
                    convolutionIndex = convolutionIndex + 1;
                end
                outputSignal(resampledIndex) = sample;
            end
            currentStep = currentStep + 1;
            resampledIndex = resampledIndex + 1;
        end
        windowIndex = windowIndex + 1;
    end

    y = outputSignal;
end