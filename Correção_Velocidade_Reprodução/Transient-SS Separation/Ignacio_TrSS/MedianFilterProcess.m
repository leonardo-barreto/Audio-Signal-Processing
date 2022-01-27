function [ x_ss,x_res,x_t] = MedianFilterProcess(senal,WINDOW,HOP,Nss,Nt,Niter)
%% function [ x_ss,x_res,x_t] = MedianFilterProcess(senal,WINDOW,NOVERLAP,Nt,Ni,Niter)
%% Input paramteres
% WINDOW window length for spectrogram computing
% HOP hop size for spectrogram computing
% Nss filter length of steady-state median filter
% Nt  filter length of transient median filter
% Niter number of iterations
%% Outputs
% x_ss steady state component signal
% x_t transient component signal
% x_res residual component signal
%%
Y=senal(:,1); %parse stereo to mono signals
[Modulo,phase]=pv_analyze(Y,WINDOW,HOP); %Spectrogram computation
[ ModuloTo,ModuloRes,ModuloIm ] = iterative_median_filter(Modulo,Nss,Nt,Niter); %Filtering stage. If the third component is 0 the separation is two component: Steady-state and Transient. Therefore ModuloRes=zeros(size(Modulo))
x_ss=pv_synthesize(ModuloTo,phase,WINDOW,HOP,HOP);%Steady-state signal resynthesize 
x_res=pv_synthesize(ModuloRes,phase,WINDOW,HOP,HOP);%Residual signal resynthesize
x_t=pv_synthesize(ModuloIm,phase,WINDOW,HOP,HOP);%Transient signal resynthesize

end

