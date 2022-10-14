function Gen_RelaxedKernel(spectrg,nFilter,time_curr)

%nFilter = 5
%time_curr = 25

[freqHeight,timeLength] = size(spectrg);


freq_c = freq;
x_c = nFilter+1;

% Kernel Parameters
s = 2; % kernel step height
l = 1; % kernel step length
h = 1 + ceil(nFilter/l)*s; % kernel total height

%kernel = zeros(size(spectrg_SS));
max_vec = zeros(2*nFilter+1,1);

%spectrg_augment = zeros(size(spectrg)+2*nFilter);
spectrg_augment = padarray(spectrg,[nFilter,nFilter],'replicate');
time_c = time_curr+nFilter;

for freq_c = 1+nFilter:freqHeight+nFilter

    for x = 0:nFilter
        h_curr = ceil(x/l)*s;
        %kernel(freq_c-h_curr:freq_c+h_curr,x_c+x) = 1;
        %kernel(freq_c-h_curr:freq_c+h_curr,x_c-x) = 1;
        max_vec(x_c+x,1) = max(spectrg(freq_c-h_curr:freq_c+h_curr,time_c+x));
        max_vec(x_c-x,1) = max(spectrg(freq_c-h_curr:freq_c+h_curr,time_c-x));
    end

end
