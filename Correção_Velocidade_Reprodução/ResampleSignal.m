%Funcao para reamostragem por fator arbitrario
%Entradas: Fator de reamostragem (resampleFactor), Sinal a ser reamostrado
%(inputSignal)

%Saidas: Sinal reamostrado: (outputSignal)

function y = ResampleSignal(samplingFrequency,resampleFactor,halfFilterCoeffs,inputSignal)
    
    %Definitions
    signalSize = length(inputSignal);
    filterCenter = 0;
    currentSampleIndex = 1;
    newFrequency = samplingFrequency*resampleFactor;
    newFrequency
    signalSize
    newPeriod = (1/newFrequency);

    %Outputs
    outputSignal = zeros (1,signalSize);
    
    
    %Code
    while filterCenter <= signalSize
        if (CheckInteger(filterCenter))
            outputSignal(currentSampleIndex) = inputSignal(currentSampleIndex);
        else
            counter = 0;
            sample = 0;
            gapRight = ceil(filterCenter) - filterCenter; %Distancias da sinc até a amostra mais próxima
            gapLeft = filterCenter - floor(filterCenter);
            while (counter < halfFilterCoeffs) %Calcula 2 termos por vez, do centro da sinc até as extremidades
                leftConvolutionIndex = floor(filterCenter) - counter;
                rightConvolutionIndex = ceil(filterCenter) + counter;
                if ( leftConvolutionIndex > 0) %Calcula lado esquerdo da sinc
                    sample = sample + inputSignal(leftConvolutionIndex)*sinc(gapLeft);
                end
                if (rightConvolutionIndex < signalSize) %Calcula lado direito da sinc
                    sample = sample + inputSignal(rightConvolutionIndex)*sinc(gapRight);
                end
                counter = counter + 1; % percorre os pontos de encontro entre a sinc e o sinal.
                gapLeft = gapLeft + 1;
                gapRight = gapRight + 1;
            end
            outputSignal(currentSampleIndex) = sample;
        end
        currentSampleIndex = currentSampleIndex + 1;
        filterCenter = filterCenter + newPeriod;
    end
    
    y = outputSignal;

end