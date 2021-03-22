% Sinthetic signal - - - - - - - - -
sig_name = 'sin_comb_pulse';
signal_name = 'sin_comb_pulse';

fs = 48000;
f0 = 600;
nHarm = 33;
t_i = .3;
t_f = .7;
dacay_flag = 0;
[x,t] = sin_pulse_gen(f0, fs, t_i, t_f, nHarm, dacay_flag);

plot_par = [.2, .8, 0, 800, 35];

x = x(round(fs*plot_par(1)) : round(fs*plot_par(2)));

%Adding white noise:
SNR = 50;
x = awgn(x, SNR);

redLines = [.1, .5];
% redLines = [];

% Spectrogram
N_w = [1024 2048 4096];

% CQT
bins_per_octave = [12 24 36];