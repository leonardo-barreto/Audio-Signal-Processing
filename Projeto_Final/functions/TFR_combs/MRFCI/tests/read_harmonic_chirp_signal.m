% Harmonic chirp signal
sig_name = 'harmonic_chirp';

amplitudes   = @(x) 1./(1:x); f0 = [200, 1000];
[x, t, f2, a2] = harmonicsource(fs_ref, dur, numpartials, f0,'f0type', ...
    'sinusoidal', 'amplitudes', amplitudes(numpartials), 'verbose', ...
    1, 'inharmonicity',0);

plot_par = [0, 1, 0, 1100, 35];
redLines = []; %test signal

% Spectrogram
N_w = [1024 2048 4096];

% CQT
bins_per_octave = [12 24 36];
