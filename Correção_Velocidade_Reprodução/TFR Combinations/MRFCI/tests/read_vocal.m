% % Combined signals - - - - -
signal_name = 'opera_fem';
NumOfSources = 2;

[x1, fs_orig] = audioread('opera_fem4.wav'); x1 = resample(x1, fs, fs_orig);
x = x1(.5*fs + 1 : 1.5*fs);

