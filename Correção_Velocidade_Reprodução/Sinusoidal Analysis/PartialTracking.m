function currentTracks = PartialTracking(inputFrame,currentTracks,DEBUG)

% Gathering frame data

currentFrame = inputFrame.currentFrame;
framePeakPowers = inputFrame.spectrumPeaks;
framePeakFreqs = inputFrame.peakFrequencies;
freqComponents = inputFrame.freqComponents;
totalFrames = inputFrame.totalFrames;

% Defining parameters

maxTracksPerFrame = 100;

freqTolerance = power(2,1/24); %quarter-tone (in Hz)
maxHysteresis = 10; % in frames.
decayThresh = 1;
minLength = 10; % in frames
maxPeakFrequency = 200000; %in Hz

% ------------------------------------ Processing EXISTING Tracks ------------------------------------------
if (~isempty(currentTracks))

    % Gathering active tracks
    activeIndexes = structArrayMatch(currentTracks,'status','active');

    %Gathering asleep tracks
    asleepIndexes = structArrayMatch(currentTracks,'status','asleep');

    % Sort in order of power (descending)
    for index = 1:length(activeIndexes)
        activeTracks(index) = currentTracks(activeIndexes(index)); 
    end

    [powerSort_Tracks,powerSort_TrackIndexes] = sortStruct(activeTracks,'currentPower',-1);

end

% ------------------------------------ Gathering all frame peak information -------------------------------------

%Discarding peaks above the max allowed frequency
[peakFreqs,freqIndexes] = sort(framePeakFreqs,'ascend');

peakFreqs(peakFreqs>maxPeakFrequency) = [];
freqIndexes(length(peakFreqs)+1:length(freqIndexes)) = [];

for index = 1:length(freqIndexes)
    peakPowers(index) = framePeakPowers(freqIndexes(index));
end

%Sorting peaks in descending order of power
[peakPowers,peakIndexes] = sort(framePeakPowers,'descend');

% -------------------------------------- Allocating new tracks if allowed-------------------------------------------


totalActiveTracks = length(activeIndexes);
totalTracks = length(currentTracks);
peakIndex = 1;

while (totalActiveTracks <= maxTracksPerFrame)
    
    currentTracks(totalTracks + 1) = createNewTrack(peakPowers(peakIndex),peakFrequencies(peakIndex),currentFrame);
    
    totalTracks = totalTracks + 1;
    totalActiveTracks = totalActiveTracks + 1;

end

    

