% % Reading audio signal:
signal_name = 'opera_fem';
NumOfSources = 2;

useTheoAlphas = 0; % For testing with controlled signals
plotTheoAlphas = 0; % For testing with controlled signals

ylim_alphas = [-6, 6];

[x1,fs_orig] = audioread('opera_fem4.wav'); x1 = resample(x1, fs, fs_orig); 
x1 = x1(2*fs + 1: 3.5*fs); 

% [x2,fs_orig] = audioread('OperaCreepy.wav'); x2 = resample(x2, fs, fs_orig); 
% x2 = x2(2*fs_orig + 1: 3*fs_orig);

% Fundamental frequency of opera_fem4.wav:---------------------------------
aux   = load('opera_fem4REF.txt'); f1 = aux(:,2);
tnovo = (2+1/fs) : 1/fs : 3.5;
f1    = interp1(aux(:,1), f1, tnovo);
%--------------------------------------------------------------------------

%x = (x1 + x2)/2; %(alterado aqui)
x = x1;