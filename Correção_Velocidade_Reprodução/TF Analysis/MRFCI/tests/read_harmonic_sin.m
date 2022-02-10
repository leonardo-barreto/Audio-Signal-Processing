% % Generating sinusoidal signals:----------------------------------------------------
signal_name = 'synth_harm_sin_noise';

NumOfSources = 2;
useTheoAlphas = 1; % For testing with controlled signals
plotTheoAlphas = 1; % For testing with controlled signals

dur = 1;

numpartials1 = 20;
amplitudes   = @(x) 1./(1:x);

y_lim_vec = [0 8000];

% First signal (synthetic): vibrato with central frequency of 500Hz
f0 = [261*4, 6];
[x1, ~, f1, a] = harmonicsource(fs, dur, numpartials1, f0,'f0type', ...
    'sinusoidal', 'amplitudes', amplitudes(numpartials1), 'verbose', ...
    1, 'inharmonicity',0);

numpartials2 = 24;
% Second signal (synthetic): vibrato with central frequency of 500Hz
f0 = [261*3, 4];
[x2, ~, f2, a] = harmonicsource(fs, dur, numpartials2, f0,'f0type', ...
    'sinusoidal', 'amplitudes', amplitudes(numpartials2), 'verbose', ...
    1, 'inharmonicity',0);

switch NumOfSources
    case 1
        x = x1;
    case 2
        x = (x1 + x2)/2;
    case 3
        x = (x1 + x2 + x3)/3;
end

%Adding white noise:
SNR = 50;

x = awgn(x,SNR);