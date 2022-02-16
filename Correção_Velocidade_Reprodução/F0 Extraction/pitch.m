function [f0,sampleStamp] = pitch(x, fs,varargin)
%pitch Estimate the fundamental frequency of audio signal
%   f0 = pitch(audioIn,fs) returns estimates of the fundamental frequency
%   over time for the audio input. Columns of the input are treated as
%   individual channels.
%
%   f0 = pitch(...,'Range',RANGE) limits the search range for the pitch
%   between the specified lower and upper band edges, inclusive. Specify
%   range in Hz as a two-element row vector of increasing values. If
%   unspecified, RANGE defaults to [50,400].
%
%   f0 = pitch(...,'WindowLength',WINDOWLENGTH) specifies the analysis
%   window length used to calculate pitch. Specify the window length in
%   samples as a positive scalar. The maximum window length is 192000. If
%   unspecified, WINDOWLENGTH defaults to round(fs*0.052).
%
%   f0 = pitch(...,'OverlapLength',OVERLAPLENGTH) specifies the number of
%   samples overlap between adjacent windows. Specify the overlap length as
%   a positive scalar smaller than the window length. If unspecified,
%   OVERLAPLENGTH defaults to round(fs*(0.042)).
%
%   f0 = pitch(...,'Method',METHOD) specifies the method used to calculate
%   the pitch. Valid inputs are:
%        PEF - Pitch Estimation Filter
%        NCF - Normalized Correlation Function
%        CEP - Cepstrum
%        LHS - Log-harmonic Summation
%        SRH - Summation of Residual Harmonics
%   If unspecified, METHOD defaults to 'NCF'.
%
%   f0 = pitch(...,'MedianFilterLength',MEDIANFILTERLENGTH) applies a
%   median filter with specified window length. The median filter is a
%   post-processing technique that operates on the estimated pitch values.
%   If unspecified, MEDIANFILTERLENGTH defaults to 1 (filter disabled).
%
%   [f0,loc] = pitch(...) returns the location associated with
%   each pitch decision. The location is the most recent sample used to
%   calculate the pitch.
%
%   EXAMPLE 1: Get the pitch contour for entire speech file
%     % Use the default settings to extract a pitch contour from a speech 
%     % file. Plot the results.
%
%       [audioIn,fs] = audioread('Counting-16-44p1-mono-15secs.wav');
%       [f0,idx] = pitch(audioIn,fs);
%
%       t = (0:length(audioIn)-1)/fs;
%       t0 = (idx - 1)/fs;
%       subplot(2,1,1); plot(t,audioIn)
%       subplot(2,1,2); plot(t0,f0)
%       xlabel('Time (s)')
%       ylabel('Pitch (Hz)')
%       ylim([50 200])
%
%
%   EXAMPLE 2: Specify nondefault parameters
%     % Use the 'PEF' method and specify an 80 ms window length with a 
%     % 10 ms hop. Limit the search range to 60-150 Hz and postprocess the
%     % pitch contour with a 3-element median filter. Plot the results.
%
%       [audioIn,fs] = audioread('SpeechDFT-16-8-mono-5secs.wav');
%
%       [f0,idx] = pitch(audioIn,fs, ...
%           'Method','PEF', ...
%           'WindowLength',round(fs*0.08), ...
%           'OverlapLength',round(fs*(0.08-0.01)), ...
%           'Range',[60,150], ...
%           'MedianFilterLength',3);
%
%       t = (0:length(audioIn)-1)/fs;
%       t0 = (idx - 1)/fs;
%       subplot(2,1,1); plot(t,audioIn)
%       subplot(2,1,2); plot(t0,f0)
%       xlabel('Time (s)')
%       ylabel('Pitch (Hz)')
%       ylim([50 200])
%
%
% See also MFCC, VOICEACTIVITYDETECTOR

% Copyright 2017-2018 The Mathworks, Inc.

%#codegen

validateRequiredInputs(x,fs)

defaults = struct( ...
    'Method',            'NCF', ...
    'Range',             cast([50,400],'like',x), ...
    'WindowLength',      cast(round(fs.*0.052),'like',x), ...
    'OverlapLength',     cast(round(fs*(0.052-0.01)),'like',x), ...
    'MedianFilterLength',cast(1,'like',x), ...
    'SampleRate',        cast(fs,'like',x), ...
    'NumChannels',       cast(size(x,2),'like',x), ...
    'SamplesPerChannel', cast(size(x,1),'like',x));

params = matlabshared.fusionutils.internal.setProperties(defaults, nargin-2, varargin{:});

validateOptionalInputs(x,fs,params)

% Determine pitch
f0 = stepMethod(x,params);

% Create sample stamps corresponding to pitch decisions
hopLength   = params.WindowLength - params.OverlapLength;
numHops     = cast(floor((size(x,1)-params.WindowLength)/hopLength),'like',x);
sampleStamp = cast(((0:numHops)*hopLength + params.WindowLength)','like',x);

% Apply median filtering
if params.MedianFilterLength ~= 1
    f0 = movmedian(f0,params.MedianFilterLength,1);
end

% Trim off zero-padded last estimate
f0 = f0(1:(numHops+1),:);
end

% -------------------------------------------------------------------------
% Validate required inputs
% -------------------------------------------------------------------------
function validateRequiredInputs(x,fs)
validateattributes(x,{'single','double'}, ...
    {'nonempty','2d','real','nonnan','finite'}, ...
    'pitch','audioIn')
validateattributes(fs,{'single','double'}, ...
    {'nonempty','positive','scalar','real','nonnan','finite'}, ...
    'pitch','fs')
end

% -------------------------------------------------------------------------
% Validate optional input
% -------------------------------------------------------------------------
function validateOptionalInputs(x,fs,userInput)
N = size(x,1);
validateattributes(userInput.Range,{'single','double'}, ...
    {'nonempty','increasing','positive','row','ncols',2,'real'}, ...
    'pitch','Range')

coder.internal.errorIf(userInput.WindowLength < 1, ...
    'audio:pitch:BadWindowLength', ...
    'WINDOWLENGTH','[1,size(x,1)]','x','round(fs*0.052)');

validateattributes(userInput.WindowLength,{'single','double'}, ...
    {'nonempty','integer','positive','scalar','real','<=',192000}, ...
    'pitch','WindowLength')
validateattributes(userInput.OverlapLength,{'single','double'}, ...
    {'nonempty','integer','scalar','real'}, ...
    'pitch','OverlapLength')
validateattributes(userInput.MedianFilterLength,{'single','double'}, ...
    {'nonempty','integer','positive','scalar','real'}, ...
    'pitch','MedianFilterLength')

coder.internal.errorIf(sum(strcmp(userInput.Method,{'NCF','PEF','CEP','LHS','SRH'}))~=1, ...
    'audio:pitch:BadMethod', ...
    'NCF','PEF','CEP','LHS','SRH');
coder.internal.errorIf(userInput.WindowLength > N, ...
    'audio:pitch:BadWindowLength', ...
    'WINDOWLENGTH','[1,size(x,1)]','x','round(fs*0.052)');
coder.internal.errorIf(userInput.OverlapLength >= userInput.WindowLength, ...
    'audio:pitch:BadOverlapLength', ...
    'OVERLAPLENGTH', 'WINDOWLENGTH');

switch userInput.Method
    case 'NCF'
        coder.internal.errorIf(fs/userInput.Range(1) >= userInput.WindowLength, ...
            'audio:pitch:BadSpecifications', ...
            'NCF','fs/RANGE(1) < WINDOWLENGTH');
        coder.internal.errorIf(fs/2<userInput.Range(2), ...
            'audio:pitch:BadSpecifications', ...
            'NCF','fs/2 >= RANGE(2)');
    case 'PEF'
        coder.internal.errorIf((userInput.Range(1)<=10) || (userInput.Range(2)>=min(4000,fs/2)), ...
            'audio:pitch:BadSpecifications', ...
            'PEF','RANGE(1) > 10 && RANGE(2) < min(4000,fs/2)');
    case 'CEP'
        coder.internal.errorIf((userInput.Range(2)>=fs/2), ...
            'audio:pitch:BadSpecifications', ...
            userInput.Method,'RANGE(2) < fs/2');
        coder.internal.errorIf(round(fs/userInput.Range(1))>2^nextpow2(2*userInput.WindowLength-1), ...
            'audio:pitch:BadSpecifications', ...
            userInput.Method,'round(fs/RANGE(1)) <= 2^nextpow2(2*WINDOWLENGTH-1)');
    case 'LHS'
        coder.internal.errorIf(((userInput.Range(2)+1)*5>=fs), ...
            'audio:pitch:BadSpecifications', ...
            userInput.Method,'(RANGE(2)+1)*5 < fs');
    case 'SRH'
        coder.internal.errorIf(((userInput.Range(2)+1)*5>=fs), ...
            'audio:pitch:BadSpecifications', ...
            userInput.Method,'(RANGE(2)+1)*5 < fs');
end
% -------------------------------------------------------------------------
end

function f0 = stepMethod(x,params)
oneCast = cast(1,'like',x);
r       = cast(size(x,1),'like',x);
c       = cast(size(x,2),'like',x);
hopLength = params.WindowLength - params.OverlapLength;

numHopsFinal = ceil((r-params.WindowLength)/hopLength) + oneCast;

% The SRH method uses a fixed-size intermediate window and hop
% length to determine the residual signal.
if isequal(params.Method,'SRH')
    N       = round(cast(0.025*params.SampleRate,'like',x));
    hopSize = round(cast(0.005*params.SampleRate,'like',x));
else
    N       = cast(params.WindowLength,'like',x);
    hopSize = cast(hopLength,'like',x);
end
numHops = ceil((r-N)/hopSize) + oneCast;

% Convert to matrix for faster processing
y = zeros(N,numHops*c,'like',x);
for channel = 1:c
    for hop = 1:numHops
        temp = x(1+hopSize*(hop-1):min(N+hopSize*(hop-1),r),channel);
        y(1:min(N,numel(temp)),hop+(channel-1)*numHops) = temp;
    end
end
% Run pitch detection algorithm
switch params.Method
    case 'SRH'
        f0 = SRH(y,params);
    case 'PEF'
        f0 = PEF(y,params);
    case 'CEP'
        f0 = CEP(y,params);
    case 'LHS'
        f0 = LHS(y,params);
    otherwise %'NCF'
        f0 = NCF(y,params);
end

% Force pitch estimate inside band edges
bE = params.Range;
f0(f0<bE(1))   = bE(1);
f0(f0>bE(end)) = bE(end);

% Reshape to multichannel
f0 = reshape(f0,numHopsFinal,c);
end

% PITCH ESTIMATION FILTER (PEF) -----------------------------------
% Gonzalez, Sira and Mike Brookes. "A Pitch Estimation Filter
% robust to high levels of noise (PEFAC)." 2011 19th European
% Signal Procesing Conference (2011): 451:455
%
function f0 = PEF(y,params)

% Read in state and parameter variables
NFFT = 2^nextpow2(2*params.WindowLength-1);
nCol = size(y,2);

logSpacedFrequency = logspace(1,log10(min(params.SampleRate/2-1,4000)),NFFT)';
linSpacedFrequency = linspace(0,params.SampleRate/2,round(NFFT/2)+1)';

wBandEdges = zeros(1,numel(params.Range),'like',y);
for i = 1:numel(params.Range)
    % Map band edges to nearest log-spaced frequency
    [~,wBandEdges(i)] = min(abs(logSpacedFrequency-params.Range(i)));
end
edge = wBandEdges;

bwTemp = (logSpacedFrequency(3:end) - logSpacedFrequency(1:end-2))/2;
bw = [bwTemp(1);bwTemp;bwTemp(end)]./NFFT;

[aFilt,numToPad] = createPitchEstimationFilter(logSpacedFrequency');

% Apply Hamming window
win = coder.const(@dsp.private.designWindow,3,size(y,1),class(y));
yw = y.*repmat(win,1,size(y,2));

% Power spectrum
Y      = fft(yw,NFFT);
Yhalf  = Y(1:(NFFT/2)+1,1:nCol);
Ypower = real(Yhalf.*conj(Yhalf));

% Interpolate onto log-frequency grid
Ylog   = interp1(linSpacedFrequency,Ypower,logSpacedFrequency);

% Weight bins by bandwidth
Ylog = Ylog.*repmat(bw,1,nCol);

% NumToPad is always scalar (indexing required for codegen)
Z   = [zeros(numToPad(1),size(Ylog,2),'like',y);Ylog];

% Cross correlation
m   = max(size(Z,1),size(aFilt,1));
mxl = min(edge(end),m - 1);
m2  = min(2^nextpow2(2*m - 1), NFFT*4);

X   = fft(Z,m2,1);
Y   = fft(aFilt,m2,1);
c1  = real(ifft(X.*repmat(conj(Y),1,size(X,2)),[],1));
R   = [c1(m2 - mxl + (1:mxl),:); c1(1:mxl+1,:)];
domain = R(edge(end)+1:end,:); % The valid domain is the second half of the correlation

% Peak-picking
locs = getCandidates(domain, edge);

f0 = logSpacedFrequency(locs);
end

% CREATE PITCH ESTIMATION FILTER ----------------------------------
function [PEFFilter,PEFNumToPad] = createPitchEstimationFilter(freq)
K     = 10;
gamma = 1.8;
num   = round(numel(freq)/2);
q     = logspace(log10(0.5),log10(K+0.5),num);
h     = 1./(gamma - cos(2*pi*q));

delta = diff([q(1),(q(1:end-1)+q(2:end))./2,q(end)]);
beta  = sum(h.*delta)/sum(delta);

PEFFilter   = (h - beta)';
PEFNumToPad = find(q<1,1,'last');
end

% NORMALIZED CORRELATION FUNCTION (NCF) ---------------------------
% B.S. Atal, "Automatic Speaker Recognition Based on Pitch
% Contours." The Journal of the Acoustical Society of America 52,
% 1687 (1972).
%
function f0 = NCF(y,params)
% Read in state and parameter variables
edge = round(params.SampleRate./fliplr(params.Range));
r    = cast(size(y,1),'like',y);

% Autocorrelation
mxl = min(edge(end),r - 1);
m2  = 2^nextpow2(2*r - 1);
c1  = real(ifft(abs(fft(y,m2,1)).^2,[],1))./sqrt(m2);
Rt  = [c1(m2 - mxl + (1:mxl),:); c1(1:mxl+1,:)];

% Energy of original signal, y
yRMS = sqrt(Rt(edge(end)+1,:));

% Clip out the lag domain of based on pitch search range
lag  = Rt( (edge(end)+1+edge(1)):end, : );

% Repmat for vectorized computation
yRMS = repmat(yRMS,size(lag,1),1);

% Normalize lag domain energy by input signal energy
lag = lag./yRMS;

% Zero-pad domain so time locs can be easily interpreted.
domain = [zeros(edge(1)-1,size(lag,2));lag];

% Peak picking
locs = getCandidates(domain,edge);

% Convert lag domain to frequency
f0 = params.SampleRate./locs;

end

% CEPSTRUM (CEP) --------------------------------------------------
% Noll, Michael A. "Cepstrum Pitch Determination". The Journal of
% the Acoustical Society of America 83, 257 (1988).
%
function f0 = CEP(y,params)
NFFT = 2^nextpow2(2*params.WindowLength-1);

% Read in state and parameter variables
edge = round(params.SampleRate./fliplr(params.Range));

% Apply Hamming window
win = coder.const(@dsp.private.designWindow,3,size(y,1),class(y));
yw = y.*repmat(win,1,size(y,2));

% Cepstral domain
domain = real(ifft(log(abs(fft(yw,NFFT)).^2)));

% Peak-picking
locs = getCandidates(domain,edge);

% Convert lag domain peaks to frequency domain
f0 = params.SampleRate./locs;
end

% LOG HARMONIC SUMMATION (LHS) ------------------------------------
% Hermes, Dik J., "Measurement of Pitch by Subharmonic Summation."
% The Journal of the Acoustical Society of America 83, 257 (1988).
%
function f0 = LHS(y,params)

% Read in state and parameter variables
edge = [ceil(params.Range(1)),floor(params.Range(end))];

% Apply Hamming window
win = coder.const(@dsp.private.designWindow,3,size(y,1),class(y));
yw = y.*repmat(win,1,size(y,2));

% Log magnitude frequency domain
S = log(abs(fft(yw,round(params.SampleRate))));

% Harmonic summation
domain = ...
    S(1:1:(edge(end)+1)  ,:) + ...
    S(1:2:(edge(end)+1)*2,:) + ...
    S(1:3:(edge(end)+1)*3,:) + ...
    S(1:4:(edge(end)+1)*4,:) + ...
    S(1:5:(edge(end)+1)*5,:);

% Peak picking
f0 = getCandidates(domain,edge);
end

% SUMMATION OF RESIDUAL HARMONICS (SRH) ---------------------------
% Drugman, Thomas and Abeer Alwan. "Joint Robust Voicing Detection
% and Pitch Estimation Based on Residual Harmonics." INTERSPEECH
% (2011).
%
function f0 = SRH(y,params)
    % Round band edges and state lpc order
    edge     = [ceil(params.Range(1)),floor(params.Range(end))];
    lpcOrder = 12;
    fs       = round(params.SampleRate);
    
    % Get LPC Residual >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    % Initialization and preallocation
    resWindowLength = cast(size(y,1),'like',y);
    numCol          = cast(size(y,2),'like',y);
    NFFT            = 2^nextpow2(2*resWindowLength-1);
    resNumHops      = numCol/params.NumChannels;
    resHopLength    = round(0.005*params.SampleRate);
    inv             = zeros(resWindowLength,numCol,'like',y);
    residual        = zeros(params.SamplesPerChannel,params.NumChannels,'like',y);
    A               = zeros(numCol,lpcOrder,'double');
    
    % Apply Hann window
    win = coder.const(@dsp.private.designWindow,2,size(y,1),class(y));
    y = y.*repmat(win,1,size(y,2));
    
    % Compute LPC
    R = cast(ifft(abs(fft(y,NFFT)).^2),'double');

    for ii = 1:numCol
        Rtemp = R(1:lpcOrder,ii);
        A(ii,:) = real(levinson(Rtemp)); % Levinson requires type double
    end
    A = cast(A,'like',y);
    
    % Compute Residual
    for ii = 1:numCol
        inv(:,ii) = real(filter(A(ii,:),1,y(:,ii)));
    end
    for jj = 1:params.NumChannels
        for kk = 1:resNumHops
            idx1 = 1+resHopLength*(kk-1);
            idx2 = min(resWindowLength+resHopLength*(kk-1),params.SamplesPerChannel);
            idx = idx1:idx2;
            residual(idx,jj) = residual(idx,jj) + inv(1:min(numel(idx),resWindowLength),kk+(jj-1)*resNumHops);
        end
    end
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    % Sum residual harmonics >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    windowLength = params.WindowLength;
    hopLength    = windowLength - params.OverlapLength;
    numHops      = ceil((size(residual,1)-windowLength)/hopLength) + 1;
    
    % Segment the residual into overlapping frames
    temp = params.NumChannels*numHops;
    sig = zeros(windowLength,temp,'like',y);
    for channel = 1:params.NumChannels
        for hop = 1:numHops
            temp = residual(1+hopLength*(hop-1):min(windowLength+hopLength*(hop-1),params.SamplesPerChannel),channel);
            sig(1:min(windowLength,numel(temp)),hop+(channel-1)*numHops) = temp;
        end
    end
    
    % Apply Blackman window
    win = coder.const(@dsp.private.designWindow,7,size(sig,1),class(y));
    sig = sig.*repmat(win,1,size(sig,2));
    
    % Magnitude frequency domain of residual
    E = abs(fft(sig,fs));
    
    % Harmonic summation
    domain = zeros(edge(end),size(E,2));
    freq   = edge(1):edge(end);
    domain(freq,:) = ...
        E(freq,:)   + ...
        E(2*freq,:) - E(round(1.5*freq),:) + ...
        E(3*freq,:) - E(round(2.5*freq),:) + ...
        E(4*freq,:) - E(round(3.5*freq),:) + ...
        E(5*freq,:) - E(round(4.5*freq),:);
    
    % Peak picking
    f0 = getCandidates(domain,edge);
    % <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
end

% GET CANDIDATES (PEAK PICKING) -----------------------------------
function locs = getCandidates(domain,edge)
numCol = size(domain,2);
locs = zeros(numCol,1,'like',domain);
lower  = edge(1);
upper  = edge(end);
assert(upper<192000);
domain = domain(lower:upper,:);
for c = 1:numCol
    [~,tempLoc] = max( domain(:,c) );
    locs(c) = tempLoc;
end
locs = lower + locs - 1;
end
