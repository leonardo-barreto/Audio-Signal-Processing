function y = PlotMeasurements(varargin)

    disp(['Number of Arguments: ' num2str(nargin)]);

    count = 1;
    while count <= nargin
        disp(['Argument ' num2str(count) ': ' inputname(count)]);
        count = count+1;
    end

    signalSize = length(cell2mat(varargin(nargin)));
    count = 1;

        measurement = ErrorMeasurements(cell2mat(varargin(count)),cell2mat(varargin(nargin)));
        numberOfMeasurements = length(measurement(:,1));
        measurements = zeros(numberOfMeasurements,signalSize,(nargin-1));

        measurements(:,:,count) = measurement;
        count = count + 1;

        while count <= (nargin-1)
            measurement = ErrorMeasurements(cell2mat(varargin(count)),cell2mat(varargin(nargin)));
            measurements(:,:,count) = measurement;
            count = count + 1;
        end
    

    indexSNR = 5;
    indexRelativeError = 4;

    x = strings(1,(nargin-1));
        underscore = strfind(inputname(1),'_');
        underscorePos = underscore(length(underscore));
        plotPoints = zeros(1,(nargin-1));

        for count = 1:(nargin-1)
            x(1,count) = extractAfter(inputname(count),underscorePos);
            plotPoints(1,count) = str2num(x(1,count));
        end

    if all(measurements(indexRelativeError,:,:) == measurements(indexRelativeError,1,:))

        RelativeErrors = zeros(1,(nargin-1));

        for count = 1:(nargin-1)
            RelativeErrors(1,count) = measurements(indexRelativeError,1,count);
        end

        figure;
        plot(plotPoints,RelativeErrors*power(10,2));
        title('Erro relativo constante! Erro relativo para cada numero de encontros da sinc');
        xlabel('Encontros com a sinc');
        ylabel('Erro relativo a amplitude (%)');

    else
        figure;
        for count = 1:(nargin-1)
            plot(measurements(indexRelativeError,:,count)*power(10,2),'DisplayName',x(1,count));
            hold on;
        end
        title('Erro relativo');
        xlabel('Amostras');
        ylabel('Erro relativo a amplitude (%)');
        hold off;
        legend;
    end

    if all(measurements(indexSNR,:,:) == measurements(indexSNR,1,:))

        SNRs = zeros(1,(nargin-1));

        for count = 1:(nargin-1)
            SNRs(1,count) = measurements(indexSNR,1,count);
        end
        figure;
        plot(plotPoints,SNRs);
        title('SNR Constante! Grafico da SNR para cada numero de encontros da sinc')
        xlabel('Encontros com a sinc');
        ylabel('SNR(dB)');
    else
        figure;
        for count = 1:(nargin-1)
           plot(measurements(indexSNR,:,count));
           hold on;
        end
        title('Razao Sinal-Ruido');
        xlabel('Numero de encontros com a sinc');
        ylabel('SNR(dB)');
        hold off;
    end
    

    y = measurements;
end