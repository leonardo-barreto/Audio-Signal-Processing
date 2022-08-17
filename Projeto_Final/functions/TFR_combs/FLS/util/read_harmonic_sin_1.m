% % Generating sinusoidal signals:----------------------------------------------------
signal_name = 'synth_harm_sin_noise';

NumOfSources = 1;
useTheoAlphas = 1; % For testing with controlled signals
plotTheoAlphas = 1; % For testing with controlled signals

dur = 1;

numpartials1 = 28;
amplitudes   = @(x) 1./(1:x);

% First signal (synthetic): vibrato with central frequency of 500Hz
f0 = [261*2, 6];
[x1, ~, f1, a] = harmonicsource(fs, dur, numpartials1, f0,'f0type', ...
    'sinusoidal', 'amplitudes', amplitudes(numpartials1), 'verbose', ...
    1, 'inharmonicity',0);

x = x1;

%Adding white noise:
SNR = 45;

x = awgn(x,SNR);