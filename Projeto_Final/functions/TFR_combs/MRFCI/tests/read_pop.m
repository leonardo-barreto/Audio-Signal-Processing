% % Reading audio signal:

signal_name = 'piano_vocal';

[x1,fs_orig] = audioread('piano_vocal.wav'); 
x1 = (x1(:,1) + x1(:,2)) / 2;
x1 = resample(x1, fs, fs_orig);

x = x1;