function [trackedFreqs] = BackwardsTracking(freqBuffer, magBuffer, param)
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

% First version - backwards peak tracking with no memory and mininum track
% length

% fs = param.fs;
% hop = param.hop;
% MSL_pd = param.MSL_pd;
% mergeNoteGap = param.mergeNoteGap;
% minNoteLength = param.minNoteLength;
% silentFrames = param.silentFrames;
nFrames = param.nFrames;
% maxNumPitches = param.maxNumPitches;
maxNumPitches = 100;
% converting the input data into a row-wise representation
freqBuffer = freqBuffer';
magBuffer = magBuffer';

% creating the output variables
trackedFreqs = zeros(nFrames,maxNumPitches);    % column track set
trackedMags = zeros(nFrames,maxNumPitches);    % column track set

%% In the first frame, all peaks should start a track
firstFrameFreqs = freqBuffer(1,:);
firstFrameMags = magBuffer(1,:);
trackedFreqs(1,1:length(firstFrameFreqs(firstFrameFreqs~=0))) = firstFrameFreqs(firstFrameFreqs~=0);
trackedMags(1,1:length(firstFrameFreqs(firstFrameFreqs~=0))) = firstFrameMags(firstFrameFreqs~=0);

%% Tracking from frame 2 to nFrames
for curFrame = 2:nFrames
    
    % selecting the current frame frequencies and magnitudes
    curFreqs = freqBuffer(curFrame,:);  % already sorted by magnitude
    curMags = magBuffer(curFrame,:);    % already sorted by magnitude
    
    % selecting the frequency indices of previous and current frames
    curFreqInd = find(curFreqs);                        % current
    pastFreqInd = find(trackedFreqs(curFrame-1,:));     % previous
    
    newTracks = zeros(1, size(trackedFreqs,2)) - 1;   % initialise new tracks with -1
    curFreqst = curFreqs(curFreqInd);               % temporary variable to allocate the current peak frequencies
    
    if ~isempty(pastFreqInd)    % if there exists available tracks in previous frame
        for i = 1:length(curFreqInd)
            if isempty(pastFreqInd)
                break;
            end
            % decision algorithm
            
        end
    end
    
end

% % List all estimated F0s and their frames
% F0Info = ListEstF0(freqBuffer, magBuffer);
% 
% % Find pairwise must-link constraints and form must-link groups
% param = {};
% param.MSL_pd = MSL_pd;                                               % Frequency deviation threshold to decide must links
% AllMustLink = FindAllMustLinks(F0Info, param);                   % Find all must-links that the domain knowledge can provide
% PreGroup = FormGroupFromMustLink(AllMustLink);                      % Form groups from must-links
% 
% % Merge groups into notes and remove small groups
% param = {};
% param.MSL_pd = MSL_pd;
% param.mergeNoteGap = round(mergeNoteGap*fs/1000/hop);                   % change to frame number
% param.minNoteLength = round(minNoteLength*fs/1000/hop);                 % change to frame number
% EstNote = FormingNotesForSomeTrack(0, PreGroup, zeros(length(PreGroup),1), F0Info, param);
% 
% % Back to F0 time-frequency representation
% param = {};
% % param.silentFrames = silentFrames;
% param.trackNum = maxNumPitches;
% param.nFrames = nFrames;
% refinedPitches = Notes2F0TimeFreqMat(EstNote, param);

