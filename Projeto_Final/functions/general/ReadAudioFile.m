function [x,signalName] = ReadAudioFile(filePath, fs_new)

    [data, fs_orig] = audioread([filePath]);
    [~,signalName,~] = fileparts(filePath);
    signalName = char(signalName);

    if size(data,2) > 1
        x = mean(data.');
    else
        x = data;
    end
    x = x(:);

    if fs_new ~= fs_orig
        x = resample(x, fs_new, fs_orig);
    end
    
end