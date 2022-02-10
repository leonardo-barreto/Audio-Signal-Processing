[data, fs_orig] = audioread([signal_name '.wav']);

if size(data,2) > 1
    x = mean(data.');
else
    x = data;
end
x = x(:);

%x = resample(x, fs, fs_orig);
%x = x(1: 15*fs);
inputSignal = x;

ylim_vec = [0 4000];
