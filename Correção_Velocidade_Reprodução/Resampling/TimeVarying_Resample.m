%Function for time-varying resampling

%Inputs: inputSignal - Signal to be resampled.
%        sinAnalysisParameters - Time-frequency analysis parameters (This should be a struct containing the fields in "Gathering Sinuosidal Analysis Data").
%        resampleFactors - Vector containing the resampling factors curve, one value per signal frame and centered around 1 (representing the original sampling rate).
%        filterCoeffs - Number of coefficients used for the sinc reconstruction filter.
%        DEBUG - Debug messages: 0 (no) or 1 (yes).

%Outputs: Resampled signal (outputSignal)

function outputSignal = TimeVarying_Resample(inputSignal,sinAnalysisParameters,resampleFactors,filterCoeffs,DEBUG)

    fprintf('\n\n------- TIME-VARYING RESAMPLING STARTED ------\n\n');

    %Gathering Sinusoidal Analysis data.

        originalSamplingRate = sinAnalysisParameters.samplingRate;              % Original sampling rate.
        timeInstants = sinAnalysisParameters.timeInstants;                      % Time instants corresponding to the center of each frame.
        frameSize = sinAnalysisParameters.windowSize;                           % Size of each frame, in samples.
        totalFrames = sinAnalysisParameters.totalFrames;                        % Total number of frames in the signal.
        hopSize = sinAnalysisParameters.hopSize;                                % Hop size between frames, in samples.

        if (length(resampleFactors) ~= totalFrames)
            error('Resampling curve size must be equal to the total number of frames, i.e. must have one factor per signal frame.\n')
        end

    %Original signal parameters
    
        inputSignalSize = length(inputSignal);
        originalPeriod = 1/originalSamplingRate; % Sampling period of the original signal.

        if DEBUG == 1
            X = sprintf('Input Signal Size: %i', inputSignalSize);
            disp(X)
            X = sprintf('Frame Size: %i', frameSize);
            disp(X)
            X = sprintf('Total Frames: %i\n\n', totalFrames);
            disp(X)
            pause(1.0);
        end

    %Preprocessing the signal for resampling: the first and last frames have to be accounted for.
        centerSamples = [1 timeInstants*originalSamplingRate inputSignalSize];

    %Outputs
        outputSignal = []; % To prevent undefined variable errors.

    %Resampling

        % Time-varying resampling is done in blocks; each frame's center sample represents the starting point of a block, and it's sampling curve factor acts on that block alone. 
        % So, usually the resampling block corresponds to the space between two center samples (roughly, the hop size).
        % The exceptions are the first/last blocks, that start/end in the first/last samples of the whole signal. Because of this, the number of blocks is one greater than the number of frames. 
        % To make transitions between resampling factors smooth, a linear interpolation is computed, raising/lowering the resampling factor slowly through the block until the next factor/block comes.
        
        totalBlocks = totalFrames+1;                % Because of edge frames, the number of blocks is the number of frames plus one.
        blockIndex = 1;                             % Index of current resampling block.
        newSamplePosition = 1;                      % New sample position: moves proportional to time, but is normalized by the original period.
        resampledOutputIndex = 1;                   % Output signal index: moves without respect to time, just the overall samples.

        while (blockIndex <= totalBlocks)

            blockSize = centerSamples(blockIndex+1)-centerSamples(blockIndex)-1;        % Blocks are disjoint, i.e. don't include the sample that starts the next one.

            % this function generates the interpolation parameters for factor transitions, at the start of each block.
            [stepSize,blockNewSamples,currentPeriod] = GenerateNormalizedParameters(resampleFactors,blockIndex,originalPeriod,blockSize,DEBUG); 

            if DEBUG == 1
                X = sprintf('Current Frame: %i of %i', blockIndex, totalFrames);
                disp(X)
                X = sprintf('Step Size: %i', stepSize);
                disp(X)
                X = sprintf('New Samples: %i', blockNewSamples);
                disp(X)
                pause(0.5);
            end 

            blockCurrentNewSample = 1;

            if newSamplePosition > inputSignalSize 
                error('New sample position exceeds the limit (greater than signal ending point).');
            end

            if (newSamplePosition < 1)
                error('New sample position exceeds the limit (less than signal starting point)');
            end

            while (blockCurrentNewSample <= (blockNewSamples+1) && newSamplePosition <= inputSignalSize)

                %This section is the traditional sinc filter convolution.
                
                convolutionIndex = 0;                                                                               %Internal index: moves from the center of sinc filter to borders to compute multiplications.
                sample = 0;                                                                                         %Current sample's value.

                if CheckInteger(newSamplePosition)                                                                  % If newSamplePosition is an integer, it coincides with an original sample and no convolution needed.
                    outputSignal(resampledOutputIndex) = inputSignal(newSamplePosition);
                else
                    leftClosestSample = floor(newSamplePosition);
                    rightClosestSample = ceil(newSamplePosition);
                    leftGap = newSamplePosition - leftClosestSample;
                    rightGap = rightClosestSample - newSamplePosition;

                    while convolutionIndex <= round(filterCoeffs/2)
                        convolutionIndex_left = leftClosestSample - convolutionIndex;                               % From center to left
                        convolutionIndex_right = rightClosestSample + convolutionIndex;                             % From center to right

                        if (convolutionIndex_left >= 1 & convolutionIndex_left < inputSignalSize)
                            sample = sample + inputSignal(convolutionIndex_left)*sinc(-leftGap - convolutionIndex); %MATLAB's sinc zeroes itself in steps of 1 already.
                        end
                        if convolutionIndex_right <= inputSignalSize
                            sample = sample + inputSignal(convolutionIndex_right)*sinc(rightGap + convolutionIndex);
                        end
                        convolutionIndex = convolutionIndex + 1;
                    end
                    outputSignal(resampledOutputIndex) = sample;
                end
                % This section moves on to the next sample.
                newSamplePosition = newSamplePosition + currentPeriod + (blockCurrentNewSample-1)*stepSize; %Interpolates between factors if needed.
                blockCurrentNewSample = blockCurrentNewSample + 1;
                resampledOutputIndex = resampledOutputIndex + 1;
            end

            blockIndex = blockIndex + 1;
        end

        fprintf('\n\n------- TIME-VARYING RESAMPLING FINISHED ------\n\n');

end