% % Combining signals - - - - -
useTheoAlphas = 0; % For testing with controlled signals
plotTheoAlphas = 0; % For testing with controlled signals

signal_name = 'violins_vibrato_drums2(4)';
NumOfSources = 3;

% [x1, fs_orig] = audioread('violin_vibrato_mix.wav'); x1 = resample(x1, fs, fs_orig);
% x1 = (x1(:,1)+x1(:,2))/2; x1 = x1(1 : 1*fs_orig);

% [x2, fs_orig2] = audioread('Drums_cut.wav'); x2 = resample(x2, fs, fs_orig2); 
% x2 = (x2(:,1)+x2(:,2))/2; x2 = x2(1 : 1*fs_orig)*.6;

% x = (x1 + x2/2)/2;
% x = x(1 : round(end/1.5));


[x1, fs_orig] = audioread('violin_drum2(4).wav'); x1 = resample(x1, fs, fs_orig);
x = (x1(:,1)+x1(:,2))/2;

