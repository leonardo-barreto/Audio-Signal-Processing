function y = PlotMeasurements(varargin)

    disp(['Number of Arguments: ' num2str(nargin)]);

    count = 1;
    while count <= nargin
        disp(['Argument ' num2str(count) ': ' inputname(count)]);
        count = count+1;
    end
   
    count = 1;
    measurements = zeros((nargin-1),nargin(ErrorMeasurements));
    while count <= (nargin-1)
        measurements(count,:) = ErrorMeasurements(varargin(count),varargin(nargin)); 
        count = count + 1;
    end

    y = measurements;

    x = zeros (1,(nargin-1));
    underscore = strfind(inputname(count),'_');
    name = 'uanisubaibf';
    underscore = underscore(length(underscore));
    for i=1:(nargin-1)
        for j=1:
        x(i) = inputname(count)(underscore)
    end
    

    
end