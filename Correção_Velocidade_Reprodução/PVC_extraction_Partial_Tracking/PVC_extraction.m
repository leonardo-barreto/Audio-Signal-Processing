% ********************************************************************* %
% Pitch Variation Curve Extraction using Partial Tracking.
% 
% Author: Ignacio Irigaray
%
% Partial Tracking Code from:
%
% J. Neri and P. Depalle, "Fast Partial Tracking of Audio with Real-Time
% Capability through Linear Programming", In Proceedings of the 21st
% International Conference on Digital Audio Effects (DAFx-18), Aveiro,
% Portugal, pp. 325-333, Sep. 2018.
%
% Ignacio Irigaray
% Universidad de la República, Montevideo, Uruguay
% ********************************************************************* %

clc;clear;close all;
addpath('utilities')

% Input audio
input_filename = 'Paulistana1.wav'; %% high pass filtered

[x, fs] = audioread(input_filename);
x = x(2*fs: 4*fs);
x=resample(x,1,4);
fs=fs/4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis Window Length
N = 2^nextpow2(.025*fs)*2-1; 
% Oversampling Factor
OverSample = 2;
% Hop Size Factor (HopSize = N/HopFactor)
HopFactor = 16;
% Magnitude Threshold for Peak-Picking (dB)
Peak_dB = -60;
% Polynomial Order Q of the short-term sinusoidal model
% 1 = Frequency, Damping
% 2 = Frequency Derivative, Damping Derivative
% 3 = Frequency 2nd Derivative, Damping 2nd Derivative, and so on..
Q = 2;
% These parameters control the assignment cost function
delta = .2; % 0<delta<1
zeta_f = 50; % Hz
zeta_a = 15; % dB
%%%%%%%%%%%%%%%%%%%%%%%%% PARTIAL TRACKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Partials, time, padding, L, S] = jun_track_partials(x,fs,...
    N,HopFactor,OverSample,Peak_dB,Q,delta,zeta_f,zeta_a);
disp('Completed tracking');
%% Can also read partials from file
[Partials, time, padding, L, fs, num_tracks] = jun_read_partials('demo_partials.bin');
%%%%%%%%%%%%%%%%%%%%%%%%% PARTIAL SYNTHESIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Pitch Variation Curve Extraction
track_F_mean=jun_ENF(Partials,time-padding(1),fs, num_tracks);
track_F_mean(track_F_mean~=track_F_mean)=1;
labels=[time/fs track_F_mean];
disp('Completed PVC extraction');
figure
plot(time/fs,track_F_mean)
title('Pitch Variation Curve')
xlabel('time (s)')
ylabel('Normalized pitch variation curve')
csvwrite('PaulisLabel.csv',labels);




