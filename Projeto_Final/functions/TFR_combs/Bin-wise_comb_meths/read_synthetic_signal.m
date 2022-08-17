% Sinthetic signal - - - - - - - - -
sig_name = 'sin_comb_pulse';
fs = 48000;
f0 = 150;
nHarm = 6;
t_i = .3;
t_f = .7;
dacay_flag = 0;
[x,t] = sin_pulse_gen(f0, fs, t_i, t_f, nHarm, dacay_flag);

%Adding white noise:
SNR = 60;
x = awgn(x,SNR);

plot_par = [.2, .8, 0, 800, 35];

% redLines = [t_i, t_f];
redLines = [];

% Spectrogram
N_w = [1024 2048 4096];

% CQT
bins_per_octave = [12 24 36];

