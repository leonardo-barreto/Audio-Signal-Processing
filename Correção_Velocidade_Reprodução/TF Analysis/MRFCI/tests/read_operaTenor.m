% % Reading audio signal:

signal_name = 'operaTenor';

NumOfSources = 1;
useTheoAlphas = 0; % For testing with controlled signals
plotTheoAlphas = 0; % For testing with controlled signals

[x1,fs_orig] = audioread('OperaTenor.wav'); x1 = resample(x1, fs, fs_orig); 
x1 = x1(2*fs + 1: 3*fs);

x = x1;