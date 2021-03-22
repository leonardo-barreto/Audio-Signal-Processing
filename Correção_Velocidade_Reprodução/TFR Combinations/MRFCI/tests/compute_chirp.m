% Sinthetic signal - - - - - - - - -
signal_name = 'chirp';

fs = 48000;
t = 0:1/fs:3.5;

t_1 = 2;
t_2 = 2.5;
t_3 = 3;

% x1 = [chirp(0:1/fs:t_1 , 100, t_1, 20000, 'logarithmic') + ...
%         chirp(0:1/fs:t_1 , 100, t_1, 20000, 'quadratic', 0, 'convex') ...
%         zeros(1, length(t_1 + 1/fs : 1/fs : t_3))];
% 
% x2 = [chirp(0:1/fs:t_2, 100, t_2, 20000, 'logarithmic') + ...
%         chirp(0:1/fs:t_2, 100, t_2, 20000, 'quadratic', 0, 'convex') ...
%         zeros(1, length(t_2 + 1/fs : 1/fs : t_3))];
% 
% x3 = chirp(0:1/fs:t_3, 100, t_3, 20000, 'logarithmic') + ...
%        chirp(0:1/fs:t_3, 100, t_3, 20000, 'quadratic', 0, 'convex');

x1 = [chirp(0:1/fs:t_1 , 100, t_1, 20000, 'logarithmic') ...
        zeros(1, length(t_1 + 1/fs : 1/fs : t_3))];

x2 = [chirp(0:1/fs:t_2, 100, t_2, 20000, 'logarithmic') ...
        zeros(1, length(t_2 + 1/fs : 1/fs : t_3))];

x3 = chirp(0:1/fs:t_3, 100, t_3, 20000, 'logarithmic');

x = x1 + x2 + x3;

% Adding white noise:
SNR = 50;
x = awgn(x,SNR);

y_lim_vec = [0 20];
