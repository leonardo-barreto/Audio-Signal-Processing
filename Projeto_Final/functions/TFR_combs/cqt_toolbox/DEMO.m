%% Demonstration
%Illustration of the use of programs to compute the
%forward and inverse constant-Q transform
%
%Christian Sch�rkhuber, Anssi Klapuri 2010-06

%% init values for CQT
fs = 44100;
bins_per_octave = 24;
fmax = fs/5;     %center frequency of the highest frequency bin 
fmin = fmax/512; %lower boundary for CQT (lowest frequency bin will be immediately above this): fmax/<power of two> 

%% generate/read input signal%
%x = randn(30*fs,1);
signal_name = '..\..\..\audio_src\paulistana3_5s.wav';

[data, fs_orig] = audioread([signal_name]);
if size(data,2) > 1
    x = mean(data.');
else
    x = data;
end
x = x(:);

if fs ~= fs_orig
    x = resample(x, fs, fs_orig);
end

% Drop frequencies outside [fmin fmax] to allow calculating 
% the SNR after inverse transform
%x = [zeros(500,1); x; zeros(500,1)];
%w1 = 2*(fmin/(fs/2)); w2 = 0.8*(fmax/(fs/2));
%[B,A] = butter(6,[w1 w2]); x = filtfilt(B,A,x); 

%% CQT
%Xcqt = cqt(x,fmin,fmax,bins_per_octave,fs);

%***computing cqt with optional input parameters***********
%Xcqt = cqt(x,fmin,fmax,bins_per_octave,fs,'q',1,'atomHopFactor',0.25,'thresh',0.0005,'win','sqrt_blackmanharris');

%***computing rasterized complex coefficients***************
 Xcqt = cqtPerfectRast(x,fmin,fmax,bins_per_octave,fs,'q',1,'atomHopFactor',0.25,'thresh',0.0005,'win','sqrt_blackmanharris');

%***precomputing filter and kernel**************************
% [B A] = butter(6,0.5,'low'); %design f_nyquist/2-lowpass filter
% K = genCQTkernel(fmax,bins_per_octave,fs,'atomHopFactor',atomHopFactor);
% Xcqt = cqt(x,fmin,fmax,bins_per_octave,fs,'kernel',K,'coeffB',B,'coeffA',A);

%***precomputing filter and kernel using cqtPerfectRast*****
% [B A] = butter(6,0.5,'low'); %design f_nyquist/2-lowpass filter
% K = genCQTkernel(fmax,bins_per_octave,fs,'atomHopFactor',atomHopFactor,'perfRast',1);
% Xcqt = cqtPerfectRast(x,fmin,fmax,bins_per_octave,fs,'kernel',K,'coeffB',B,'coeffA',A);

%% inverse CQT
%y = icqt(Xcqt);
%SNR = 10*log10(sum(abs(x).^2)/sum(abs(x-y).^2)); %calculate the SNR

%% plot CQT
%plotCQT(Xcqt,fs,1,'surf');

%% get CQT (interpolated)
%intCQT = getCQT(Xcqt,'all','all');
absCQT = abs(Xcqt.spCQT);


 plot_range = 100;
 plot_max = max(max(10*log10(absCQT)));
 emptyHops = Xcqt.intParams.firstcenter/Xcqt.intParams.atomHOP;
 maxDrop = emptyHops*2^(Xcqt.octaveNr-1)-emptyHops;
 droppedSamples = (maxDrop-1)*Xcqt.intParams.atomHOP + Xcqt.intParams.firstcenter;
 t = (1:size(absCQT,2))*Xcqt.intParams.atomHOP-Xcqt.intParams.preZeros+droppedSamples;
 t = t./fs;
 f = 1:size(absCQT,1);
 PlotSpectrogram_ylin(f,t,[plot_max-plot_range plot_max],10*log10(absCQT));
 title(sprintf('Espectrograma CQT'))
 set(gca,'YTick',1:Xcqt.bins/2:Xcqt.octaveNr*Xcqt.bins);
 h = get(gca); yTick = h.YTick';
 yTickLabel = num2str(round(Xcqt.fmin*2.^((yTick-1)/Xcqt.bins)),5);
 set(gca,'YTickLabel',yTickLabel);


