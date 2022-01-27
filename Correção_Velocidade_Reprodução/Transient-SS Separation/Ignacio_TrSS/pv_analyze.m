function [Moduli,Phases]=pv_analyze(X,win,hop)

[nr, nc] = size (X);
% compute the window coefficients
%
WIN_COEF = hanningz(win);
%WIN_COEF = hanning(win);

%
% create a matrix Z whose columns contain the
% windowed time-slices
%
num_win =ceil((nr-win+hop)/hop);
%Z= zeros (win,num_win);
%
% now modulate the signal with the window,
% frame by frame taking due care that the
% signal is zero padded at the end
%
start = 1;
textprogressbar('Analiyse: ')

for i = 0:num_win
    frame_end = win-1;
    if (start+frame_end >= size(X,1))
         frame_end = size(X,1)-start;
    end
    win_end = frame_end+1;
    Z = X(start:start+frame_end).*WIN_COEF(1:win_end);
    FZ(1:win_end,i+1) = fft(fftshift(Z));
    start = start + hop;
textprogressbar(i/num_win*100);
end; 

Moduli = abs(FZ);
Phases = angle(FZ);
textprogressbar('Done')

end