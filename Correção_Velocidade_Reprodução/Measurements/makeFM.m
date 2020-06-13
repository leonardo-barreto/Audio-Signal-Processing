function [outputSignal,modulationCurve] = makeFM(k,duration)


omega0 = 2000; %rad/s
omegaV = 10; %rad/s
fs = 40000; % in Hz
alfa = 0.01;

n = 0:1:(duration*fs)-1;
outputSignal = [];

modulationCurve = 1 + alfa*cos(omegaV*n./fs);

outputSignal = cos(omega0*k*(n./fs + (alfa/omegaV)*sin(omegaV*n./fs)))/k^5;

plot(n,modulationCurve,'LineWidth',2);
X = sprintf('Curva de modulação');
title(X,'FontSize', 18);
xlabel('Tempo (s)','FontSize', 18);
ylabel('Modulação da portadora (Hz)','FontSize', 18);



end