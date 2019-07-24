%Funcao para reamostragem por fator arbitrario
%Entradas: Fator de reamostragem (resampleFactor), Sinal a ser reamostrado
%(inputSignal)

%Saidas: Sinal reamostrado: (outputSignal)

function y = ResampleSignal(inputSignal,resampleFactor,filterCoeffs)
    
    %Definitions
    inputSignalSize = length(inputSignal);
    newSamplePosition = 1;
    resampledIndex = 1;
    newPeriod = (1/resampleFactor);

    %Outputs
    outputSignal = zeros(1,inputSignalSize); % Para prevenir erros de undefined variable
   
    %Code
    while newSamplePosition <= inputSignalSize % percorre o sinal de entrada no novo período de amostragem
            convolutionIndex = 1; % índice interno que percorre o sinal de entrada para a multiplicação pela sinc
            sample = 0;
            if CheckInteger(newSamplePosition) % se a amostra na nova taxa coincide com a amostra na taxa antiga
                outputSignal(resampledIndex) = inputSignal(newSamplePosition);
            else
                
                leftGap = newSamplePosition - floor(newSamplePosition);
                rightGap = ceil(newSamplePosition) - newSamplePosition;

                while convolutionIndex <= round(filterCoeffs/2)
                    convolutionIndex_left = (newSamplePosition - leftGap) - (convolutionIndex-1);
                    convolutionIndex_right = (newSamplePosition + rightGap) + (convolutionIndex-1);
                    if convolutionIndex_left > 0
                        sample = sample + inputSignal(convolutionIndex_left)*sinc(-leftGap - (convolutionIndex-1));
                    end
                    if  convolutionIndex_right < inputSignalSize
                        sample = sample + inputSignal(convolutionIndex_right)*sinc(rightGap + (convolutionIndex-1));
                    end
                    convolutionIndex = convolutionIndex + 1;
                end
                outputSignal(resampledIndex) = sample;
            end
            newSamplePosition = newSamplePosition + newPeriod;
            resampledIndex = resampledIndex + 1;
    end
    y = outputSignal;
end