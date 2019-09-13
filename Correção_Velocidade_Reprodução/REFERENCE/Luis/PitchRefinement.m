function [refinedPitches, EstNote] = PitchRefinement(estPitches, magBuffer, param)
% Calculate STFT
% Input:
%   - wavData           : wave form mono
%   - para
%       - frameLen      : frame length in points
%       - hop           : hop size in points
%       - window        : a window function vector of frameLength
%       - fftLen        : fft length in points
% Output:
%   - SpecData          : complex spectrogram, dimension = [fftLength/2+1, total frame number]
%
% Author: Zhiyao Duan
% Last modified: 5/24/2013

fs = param.fs;
hop = param.hop;
MSL_pd = param.MSL_pd;
mergeNoteGap = param.mergeNoteGap;
minNoteLength = param.minNoteLength;
% silentFrames = param.silentFrames;
nFrames = param.nFrames;
maxNumPitches = param.maxNumPitches;

% List all estimated F0s and their frames
F0Info = ListEstF0(estPitches, magBuffer);

% Find pairwise must-link constraints and form must-link groups
param = {};
param.MSL_pd = MSL_pd;                                               % Frequency deviation threshold to decide must links
AllMustLink = FindAllMustLinks(F0Info, param);                   % Find all must-links that the domain knowledge can provide
PreGroup = FormGroupFromMustLink(AllMustLink);                      % Form groups from must-links

% Merge groups into notes and remove small groups
param = {};
param.MSL_pd = MSL_pd;
param.mergeNoteGap = round(mergeNoteGap*fs/1000/hop);                   % change to frame number
param.minNoteLength = round(minNoteLength*fs/1000/hop);                 % change to frame number
EstNote = FormingNotesForSomeTrack(0, PreGroup, zeros(length(PreGroup),1), F0Info, param);

% Back to F0 time-frequency representation
param = {};
% param.silentFrames = silentFrames;
param.trackNum = maxNumPitches;
param.nFrames = nFrames;
refinedPitches = Notes2F0TimeFreqMat(EstNote, param);

