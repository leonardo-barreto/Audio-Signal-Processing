close all
filen = '1609.aif'; 
[x,fs]=audioread(filen); 
x=x(:,1); %Left channel only.
filtro_enf_narrow= designfilt('bandpassfir', 'StopbandFrequency1', 88, 'PassbandFrequency1', 90, 'PassbandFrequency2', 98, 'StopbandFrequency2', 101, 'StopbandAttenuation1', 23, 'PassbandRipple', 1, 'StopbandAttenuation2', 40, 'SampleRate', fs,'DesignMethod','kaiserwin');% diseño del filtro

%% 
x_enf=filtro_enf_narrow.filtfilt(x);%filter the signal
z_enf=hilbert(x_enf);  % calculates the analytic signal
instfrq_enf = fs/(2*pi)*diff(unwrap(angle(z_enf)));% Instantaneous frequency
plot(instfrq_enf,'r')
%% 
t=(0:(length(instfrq_enf)-1))'/fs;
audiowrite('1609_enf.wav',x_enf,fs);
%% downsample ENF
f0_enf=resample(instfrq_enf,1,100);
t_enf=t(1:100:end);
dlmwrite('ENF_variation.csv',[t_enf,f0_enf]); % write enf labels



