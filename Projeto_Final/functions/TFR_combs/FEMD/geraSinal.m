function [x,ft] = geraSinal(f0,fs,dur,numpartials,pervibrato)
%GERASINAL Gera um sinal sintético harmônico x
% Parâmetros de entrada:
% f0 = [f0Inicial f0Final f0Vibrato]
% fs: frequência de amostragem
% dur = [durGlissando durVibrato]
% numpartials: número de parciais
% pervibrato: percentual de variação frequencial no vibrato

% Gerando o chirp de terceira ordem:---------------------------------------
a = 4*(f0(2)-f0(1))/dur(1)^3;
b = -6*(f0(2)-f0(1))/dur(1)^2; 
c = 3*(f0(2)-f0(1))/dur(1);
t1  = 0:1/fs:dur(1);
ft1 = @(x) a*x.^3 + b*x.^2 + c*x + f0(1);
% Gerando o chirp senoidal:------------------------------------------------
t2 = 0:1/fs:dur(2);
ft2 = sin(2*pi*f0(3)*t2); 
aux = length(0:1/fs:(dur(2)/2.5));
h = [hanning(aux); zeros(length(t2)-aux,1)] + ...
    [zeros(floor(aux/2),1); hanning(aux); zeros(aux-1,1)] + ...
    [zeros(aux-1,1); hanning(aux); zeros(floor(aux/2),1)] + ...
    [zeros(length(t2)-aux,1); hanning(aux)];
ft2 = (ft2(:).*h(:)*pervibrato/100 + 1)*ft1(dur(1)/2);
% Concatenando para obter chirp de terceira ordem + senoidal:--------------
ft  = [ft1(t1(1:floor(length(t1)/2))) ft2' ft1(t1(ceil(length(t1)/2):end))];
phi = cumsum(ft);
% Sintetizando o sinal de áudio:-------------------------------------------
amplitudes = @(x) 1./((1:x).^2); 
amp = amplitudes(numpartials)/sum(amplitudes(numpartials)); 
x   = 0;
for n=1:numpartials
    x = x + amp(n)*sin(2*pi*n*phi/fs);
end
x = x(:);
end

