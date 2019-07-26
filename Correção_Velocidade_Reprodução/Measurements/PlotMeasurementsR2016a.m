function y = PlotMeasurements(varargin)

    disp(['Number of Arguments: ' num2str(nargin)]);

    count = 1;
    while count <= nargin
        disp(['Argument ' num2str(count) ': ' inputname(count)]);
        count = count+1;
    end

    signalSize = min (length(cell2mat(varargin(1))),length(cell2mat(varargin(nargin))));

    disp(['Signal Size: ' num2str(signalSize)]);
    disp('Making Measurements...');
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

    disp('Measurements done!');
    disp(['Signals evaluated: ' num2str(length(measurements(1,1,:)))]);
    disp(['Measurements evaluated: ' num2str(length(measurements(:,1,1)))]);
    disp(['Signal Size: ' num2str(length(measurements(1,:,1)))]);

    indexSNR = 5;
    indexRelativeError = 4;


    plotPoints = [128,256,512,1024];
    x = ['0128';'0256';'0512';'1024'];

    numberOfSignals = (nargin-1);
    constant = 1;
    signals = 1;

    while constant == 1 & signals <= numberOfSignals
        if ~all(measurements(indexRelativeError,:,signals))
          constant = 0;
        end
        signals = signals + 1;
    end


    if constant == 1

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
            plot(measurements(indexRelativeError,:,count)*power(10,2),'DisplayName',x(count,:));
            hold on;
        end
        title('Erro relativo');
        xlabel('Amostras');
        ylabel('Erro relativo a amplitude (%)');
        hold off;
        legend;
    end

    constant = 1;
    signals = 1;

    while (constant == 1 & signals <= (nargin-1))
        if ~all(measurements(indexSNR,:,signals))
          constant = 0;
        end
      signals = signals + 1;
    end

    if constant == 1

        SNRs = zeros(1,(nargin-1));

        for count = 1:(nargin-1)
            SNRs(1,count) = measurements(indexSNR,1,count);
        end
        figure;
        plot(plotPoints,SNRs);
        title('SNR Constante! SNR para cada numero de encontros da sinc')
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
