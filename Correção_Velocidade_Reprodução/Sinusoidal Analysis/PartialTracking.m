function y = PartialTracking(inputFrame,currentTracks,templateTrack,MAXTRACKS,DEBUG)

% STARTING TRACK PARAMETER CANDIDATES
%
% FREQUENCY TOLERANCE
% STARTING AMPLITUDE
% 
%
%

% ENDING TRACK PARAMETER CANDIDATES
%
% MAGNITUDE DECAY
% LENGTH OF TRACK
% CURRENT AMPLITUDE
% 
%
%

% Gathering frame data
    
framePeaks = inputFrame.spectrumPeaks;
framePeakPositions = inputFrame.peakPositions;
freqComponents = inputFrame.freqComponents;
totalFreqBins = inputFrame.totalFreqBins;
totalFrames = inputFrame.totalFrames;
currentFrame = inputFrame.currentFrame;
totalTracks = inputFrame.totalTracks;

% Defining parameters

freqTolerance = 1;
decayThreshold = 1;
lengthThreshold = 1;

% ------------------------------------ Processing EXISTING Tracks ------------------------------------------

% Gathering active tracks
activeIndexes = StructArrayMatch(currentTracks,'status','active');

%Gathering asleep tracks
asleepIndexes = StructArrayMatch(currentTracks,'status','asleep');

% Sort in order of power (descending)
for index = 1:length(activeIndexes)
    powerSort_Tracks(index) = currentTracks(activeIndexes(index)); 
end

[powerSort_Tracks,powerSort_TrackIndexes] = sortStruct(currentTracks,'currentPower',-1);

% ------------------------------------ Gathering all peak information -------------------------------------

[powerSort_Peaks,powerSort_PeakIndexes]

% ------------------------------------ New track candidates -----------------------------------------------

while (totalTracks <= MAXTRACKS)
    

