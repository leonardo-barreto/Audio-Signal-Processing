% % Combined signals - - - - -
useTheoAlphas = 0; % For testing with controlled signals
plotTheoAlphas = 0; % For testing with controlled signals

signal_name = 'violins_vibrato';
NumOfSources = 2;

[x1, fs_orig] = audioread('violin_vibrato_mix.wav'); x1 = resample(x1, fs, fs_orig);
x1 = (x1(:,1)+x1(:,2))/2; x = x1(1 : 1*fs);


