% Set parameters for music
% "%--" before a parameter indicates that this parameter has been initialized before but used here as well
% "%--" before a comment of a parameter indicates that this parameter is important
%
% Author: Zhiyao Duan
% Created: 6/20/2012
% Last modified: 5/24/2013

% parameters for multiple pitch estimation
frameLen =8*512;                            %-- 92ms for fs=44100
fftLen = 16*1024;                             % FFT number of points (zero padding)
win = hann(frameLen, 'periodic');           % window function
hop = 256;                                  %-- 10ms for fs=44100

weightCoeff = 0.9;
bayThreshold = 0.01;
MSL_pd = 0.01;                               %-- pitch difference threshold of must-link (midi number
tonalnessTh = 0.8;                          % tonalness threshold
specThreshFactor = 0.004;                   % spectral magnitude threshold factor
mergeNoteGap = 50;                         %-- threshold for note formation (ms), two notelets with gap less than this threshold can be merged
minNoteLength = 2*50;                        %-- threshold for note length (ms), notelets shorter than this will be removed


RMS_th = sqrt(0.01);                        % frame whose RMS lower than this threshold will be considered silent
F0max = 20000;                              % maximum frequency up to which the spectrum is evaluated
F0min = 55;                                 % the range of detectable pitches starts in F0min
% alpha_at = 0.10;                             % smoothing coeff for the amplitude threshold
alpha_at = 1500/fftLen;
maxNumPitches = 100;                         % maximum number of detectable pitches per frame
specSalienceFactor = (0.1)^(0.25);          % spectral salience threshold factor
compressionParam = 0.5;                     % IFFT compression parameter
autocorrThreshFactor = 0.1;               % autocorr threshold factor
autocorrSalienceFactor = 0.25;               % autocorr salience threshold factor
finalSalienceThresh = 4e-6;                 % final salience threshold
doNoteTracking = 1;                          % 0: not refine F0 estimates again, 1: refine F0 estimates again, by removing some outliers and filling some gaps

eps_at = 0.5298;
eps_pk = 0.1443;

