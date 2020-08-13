[d,sr] = audioread('audiochirp1.wav');
S = spectrogram(d,hann(8192),6144,8192);          % SP Toolbox routine (or use ifgram.m below)
[R,M]=extractrax(abs(S));     % find peaks in STFT *magnitude*
disp(['size of R is ',num2str(size(R,1)),' rows x ',num2str(size(R,2)),' cols']);

           % default specgram step is NFFT/2 i.e. 128
F = R*sr/8192;                  % Convert R from bins to Hz
[s,f,t] = spectrogram(d,hann(8192),6144,8192,sr);

colormap(1-gray)              % black is intense, white is quiet
hold on
plot(t,F,'r');              % the tracks follow the specgram peaks