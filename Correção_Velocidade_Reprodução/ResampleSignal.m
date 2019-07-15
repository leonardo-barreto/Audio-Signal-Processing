%Funcao para reamostragem por fator arbitrario
%Entradas: Fator de reamostragem (resampleFactor), Sinal a ser reamostrado
%(inputSignal)

%Saidas: Sinal reamostrado: (outputSignal)

function y = ResampleSignal(x1,x2,x3)

    %Input Arguments
    resampleFactor = x1; %REAL
    halfFilterCoeffs = x2; %INTEGER
    inputSignal = x3; %VECTOR
    
    %Definitions
    signalSize = length(inputSignal);
    filterCenter = 0;
    currentSampleIndex = 0;
    newPeriod = 1/resampleFactor;

    %Outputs
    outputSignal = zeros (1,signalSize);
    
    %Code
    while filterCenter < signalSize
        if (CheckInteger(filterCenter))
            outputSignal(currentSampleIndex) = inputSignal(currentSampleIndex);
        else
            counter = 0;
            sample = 0;
            rightStep = ceil(filterCenter) - filterCenter;
            leftStep = filterCenter - floor(filterCenter);
            while (counter < halfFilterCoeffs)
                endSignal = 0;
                if ((index = floor(filterCenter) - counter) > 0) %Calcula lado esquerdo da sinc
                    sample = sample + inputSignal(index)*sinc();
                end

            end

        end
    

        
    end
    
    
    

end