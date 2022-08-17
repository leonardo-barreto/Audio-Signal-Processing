% % Combining signals - - - - -
signal_name = 'violins_piano';
NumOfSources = 3;

[x1, fs_orig] = audioread('violin_vibrato_mix.wav'); x1 = resample(x1, fs, fs_orig);
x1 = (x1(:,1)+x1(:,2))/2;
x1 = x1(1 : 1*fs);


[x2, fs_orig] = audioread('sonata_MIDI.wav'); x2 = resample(x2, fs, fs_orig);
x2 = x2(1 : 1*fs);

x = (x1 + x2)/2;