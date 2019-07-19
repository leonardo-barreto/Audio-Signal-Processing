%Funcao para reamostragem por fator arbitrario
%Entradas: Fator de reamostragem (resampleFactor), Sinal a ser reamostrado
%(inputSignal)

%Saidas: Sinal reamostrado: (outputSignal)

function outputSignal = ResampleSignal(resampleFactor,filterCoeffs,inputSignal)
    
    %Definitions
    signalSize = length(inputSignal);
    newIndex = 0;
    currentSampleIndex = 1;
    newPeriod = (1/resampleFactor);

    %Outputs
    outputSignal=zeros(1,signalSize);
   
    %Code
    while newIndex < signalSize
            originalIndex = 1;
            sample = 0;
            while (originalIndex <= filterCoeffs)
                sample = sample + inputSignal(originalIndex)*sinc(newIndex - (originalIndex-1));
                originalIndex = originalIndex + 1;
            end
        newIndex = newIndex + newPeriod;
        outputSignal(currentSampleIndex) = sample;
    end


end