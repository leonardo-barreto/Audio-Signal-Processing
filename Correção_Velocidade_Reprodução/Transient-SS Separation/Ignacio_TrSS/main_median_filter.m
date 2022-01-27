%%Main file for transient and steady-stae separation.
% WINDOW: window length for spectrogram computing
% HOP: hop size for spectrogram computing
% Nss: filter length of steady-state median filter
% Nt:  filter length of transient median filter
% Niter number of iterations: niter= 0 two transient and steady-state
%   decomposition, niter>0 transient, residual and steady-state decomposition
clear variables
close all
file='prueba.wav';
[A,FS,NBITS]=wavread(file);%read file
WINDOW=2048;
HOP=512;
Nss=8;
Nt=8;
Niter=1;% 
[x_ss,x_res,x_t]=MedianFilterProcess(A,WINDOW,HOP,Nss,Nt,Niter);
wavwrite(x_ss,FS,[file '_SteadyState.wav']);
wavwrite(x_t,FS,[file '_Transient.wav']);
wavwrite(x_res,FS,[file '_Residual.wav']);