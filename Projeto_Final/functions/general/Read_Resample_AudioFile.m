function x = Read_Resample_AudioFile(signal_name, fs)

    [data, fs_orig] = audioread([signal_name]);

    if size(data,2) > 1
        x = mean(data.');
    else
        x = data;
    end
    x = x(:);

    if fs ~= fs_orig
        x = resample(x, fs, fs_orig);
    end

    fs_out = fs;

end