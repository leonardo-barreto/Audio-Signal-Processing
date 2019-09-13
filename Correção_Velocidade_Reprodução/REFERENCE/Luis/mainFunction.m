clc; clear;
% close all;
% main function

%% Set parameters
SetParameters;
tic
%% Section I - Loading the wav file
wavPath = 'sounds\';                        % path where the sound file is located
wavFile = 'orchestraSmallSin';                          % sound file name
fileFullPath = [wavPath wavFile '.wav'];    % full path to access the file
[wavData, fs] = audioread(fileFullPath);    % load the file
wavData = mean(wavData, 2);                 % convert stereo to mono (optional)
wavData = wavData - mean(wavData);          % zero mean (optional)

%% Section II - Obtaining the spectral representation

% Short-time Fourier transform (with normalisation in case of deterministic peak detection) 
param = {};
param.frameLen = frameLen;
param.hop = hop;
param.win = win;
param.fftLen = fftLen;

[magSpecData, nFrames] = CalculateSTFT(wavData, param, 0);

% Plot spectrogram (optional)
param = {};
param.fs = fs;
param.nFrames = nFrames;
param.fftLen = fftLen;

PlotSpectrogram (magSpecData, param);

%% Section III - Peak detection (deterministic and bayesian cases)

% Deterministic: tonalness-based
param = {};
param.frameLen = frameLen;
param.fftLen = fftLen;
param.fs = fs;
param.nFrames = nFrames;
param.alpha_at = alpha_at;
param.specThreshFactor = specThreshFactor;
param.tonalnessTh = tonalnessTh;

freqBufferInput = zeros(maxNumPitches,nFrames);

[freqBuffer, magBuffer] = TonalnessPeakDetection (magSpecData, freqBufferInput, param);

% Bayesian peak detection:
param = {};
param.frameLen = frameLen;
param.fftLen = fftLen;
param.fs = fs;
param.hop = hop;
param.win = win;
param.nFrames = nFrames;
param.bayThreshold = bayThreshold;

% [freqBuffer, magBuffer] = BayesianPeakDetection (magSpecData, wavData, freqBufferInput, param);


% Plotting spectrogram with the detected frequencies
param = {};
param.fs = fs;
param.nFrames = nFrames;
param.fftLen = fftLen;

PlotSpecWithPeaks(magSpecData, freqBuffer, param);

%% Section IV - Peak Tracking
param = {};
param.fs = fs;
param.hop = hop;
param.MSL_pd = MSL_pd;
param.mergeNoteGap = mergeNoteGap;
param.minNoteLength = minNoteLength;
param.nFrames = nFrames;
% param.silentFrames = silentFrames;
param.maxNumPitches = maxNumPitches;

[trackedFreqs, EstNote] = PitchRefinement(freqBuffer, magBuffer, param);
% [trackedFreqs, EstNote] = BackwardsTracking(freqBuffer, magBuffer, param);


% Plotting spectrogram with the tracked frequencies
param = {};
param.fs = fs;
param.nFrames = nFrames;
param.fftLen = fftLen;

PlotSpecWithPeaks(magSpecData, trackedFreqs, param);
toc
%% global pitch variation curve
[aCurve, wCurve] = GPVCurve(EstNote,nFrames, weightCoeff);
smoothCurve = movmean(wCurve+1,10);
figure; plot(aCurve);
figure; plot(smoothCurve);ylim([0.98 1.02]);
centralSamples = fix((0:nFrames-1)*hop + (frameLen+1)/2);

%% non-uniform resampling

xRestored = TimeVaryingResample(wavData, smoothCurve, centralSamples');
