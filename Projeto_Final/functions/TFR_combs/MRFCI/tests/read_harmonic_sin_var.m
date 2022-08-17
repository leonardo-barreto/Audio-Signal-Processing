% % Generating sinusoidal signals:----------------------------------------------------
signal_name = 'synth_harm_sin_var_noise';

NumOfSources = 1;

dur = 1.5;

numpartials1 = 20;
amplitudes   = @(x) 1./(1:x);

% First signal (synthetic): vibrato with central frequency of 500Hz
 f0 = [880, 4];
[x1, ~, f1, a] = harmonicsource(fs, dur, numpartials1, f0,'f0type', ...
'cresc_vib', 'amplitudes', amplitudes(numpartials1), 'verbose', ...
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

x = x';