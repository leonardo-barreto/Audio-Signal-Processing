function [max_matrix] = Gen_RelaxedKernel(spectrg,nFilter,time_curr)

    %TESTING PURPOSES
    %nFilter = 5;
    %time_curr = 25;
    %spectrg = zeros(50,50);
    %spectrg(11:20,20) = 1:10;
    %spectrg(21:32,20) = 11:22;
    %spectrg(11:20,23) = 1:10;
    %spectrg(21:32,23) = 11:22;

    % Kernel Parameters
    s = 2; % kernel step height
    l = 8; % kernel step length
    h = ceil(nFilter/l)*s; % kernel max height

    [freqHeight,timeLength] = size(spectrg);
    max_matrix = zeros(freqHeight,2*nFilter+1); % Output matrix

    if (time_curr <= nFilter | time_curr > timeLength-nFilter)
        spectrg_augment = padarray(spectrg,[h,nFilter],'replicate'); % Padding in the spectrogram to accomodate border kernels
        pad_horiz = nFilter;
    else
        spectrg_augment = padarray(spectrg,[h,0],'replicate'); % Padding in the spectrogram to accomodate border kernels
        pad_horiz = 0;
    end

    time_c = time_curr+pad_horiz;     % Kernel horizontal center (in spectrogram)
    x_c = nFilter+1;                  % Kernel horizontal center (in output vector)

    for freq_bin = 1+h:freqHeight+h
        y = freq_bin-h;
        for x = 0:nFilter
            h_curr = ceil(x/l)*s;
            %spectrg_augment(freq_bin-h_curr:freq_bin+h_curr,time_c+x) = 1; % Testing purposes
            %spectrg_augment(freq_bin-h_curr:freq_bin+h_curr,time_c-x) = 1; % Testing purposes
            max_matrix(y,x_c+x) = max(spectrg_augment(freq_bin-h_curr:freq_bin+h_curr,time_c+x));
            max_matrix(y,x_c-x) = max(spectrg_augment(freq_bin-h_curr:freq_bin+h_curr,time_c-x));
        end

    end

end
