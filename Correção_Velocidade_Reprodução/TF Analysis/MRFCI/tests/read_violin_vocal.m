% % Combined signals - - - - -
useTheoAlphas = 0; % For testing with controlled signals
plotTheoAlphas = 0; % For testing with controlled signals

signal_name = 'opera_fem';
NumOfSources = 2;

[x1, fs_orig] = audioread('opera_fem4.wav'); x1 = resample(x1, fs, fs_orig);
x = x1(.5*fs + 1 : 1.2*fs);

% Fundamental frequency of opera_fem4.wav:---------------------------------
aux = load('opera_fem4REF.txt'); f1 = aux(:,2);
tnovo = (.5 + 1/fs) : 1/fs : 1.5;
f1= interp1(aux(:,1), f1, tnovo);
%--------------------------------------------------------------------------