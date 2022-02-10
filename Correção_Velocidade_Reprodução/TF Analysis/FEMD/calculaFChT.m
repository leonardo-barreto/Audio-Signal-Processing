function [specs,mag_specs] = calculaFChT(xn, fs, N, Nalfas, alfas)
% CALCULAFCHT retorna inst?ncias da FChT para os diferentes valores de alfa
% pr?-determinados
% Sinais de entrada:
% xn: sinal de entrada
% fs: frequ?ncia de amostragem
% Nalfas: n?mero de valores de alfa
% alfas: valores de alfa 
% Os valores de N e M s?o explicados no artigo: 
% "Fan-chirp transform for music representation"
% specs: FChT para todos os bins de frequ?ncia (para cada warp)
% mag_specs: m?dulo da FChT para todos os bins (para cada warp)

% warping functions
psi_wrap = @(x,alfa) (-1+sqrt(1+2*alfa*x))/alfa;
phi_wrap = @(x,alfa) (1+1/2*alfa*x)*x;

max_alfa = max(alfas);
M = 2^nextpow2(ceil(N / (1-abs(max_alfa/fs) * N/2)));

% Nalfas espectros (todas as frequencias):
specs = zeros(M,Nalfas);
mag_specs = zeros(M,Nalfas);

for i = 1:Nalfas
    % Valor atual de alfa:
    alfa_atual = alfas(i);
    % Reamostrando:
    if alfa_atual ~= 0 
        tn = ((0:1:N-1)-(N-1)/2)/fs;
        tr = phi_wrap(tn(1),alfa_atual) + ((0:1:M-1)+1/2)*(tn(end)-tn(1))/M;
        tt = psi_wrap(tr,alfa_atual);
        xr = interp1(tn,xn,tt);
    else
        tn = ((0:1:N-1)-(N-1)/2)/fs;
        tt = tn(1) + ((0:1:M-1)+1/2)*(tn(end)-tn(1))/M;
        xr = interp1(tn,xn,tt);
    end
    % A multiplicacao pela janela eh realizada apos a deformacao no tempo:
    Xr = fft(xr'.*hanning(M));
    specs(:,i)     = Xr;
    mag_specs(:,i) = abs(Xr);
end

end