[y, fs] = audioread('b.wav');
%sound(y,fs)

% 
% title("N = 128")
% 
% subplot(2,1,2)
% s = spectrogram(y, hamming(8192), 4096, 1024, fs, 'yaxis');
%title("N = 128")

%produz o espectrograma do sinal analisado
[s,f,t,p] = spectrogram(y, hamming(512), 64, 1024, fs, 'yaxis');
spectrogram(y, hamming(128), 64, 1024, fs, 'yaxis');

%construção da matriz auxiliar de potencias
%iremos definir o melhor método de estimação de limiar 
%para as matrizes

timeLength = length(t);
freqLength = length(f);
for i = 1:timeLength
    for k = 1:freqLength
        powerFrames(i,k)=p(k,i);
    end
end

powerFramesdB = 10*log10(powerFrames);
%plot(powerFramesdB(1:timeLength,1:freqLength)') // Plot linhas e colunas
%arbitrárias

%Geração de polinomio de limiar para um único frame

%plot(powerFramesdB(1,1:freqLength)') 
w = powerFramesdB(1,1:freqLength);
polyCoef = polyfit(1:freqLength,w,3); % gera os coeficientes do polinomio que aproxima x com ordem n
fitPoly = polyval(polyCoef,1:freqLength); % gera o polinomio com os coeficientes anteriores

%Visualização de um frame com seu polinomio de limiar

plot(fitPoly)
hold on
plot(w)

%Detecção de pico de um único frame

plot(w)
hold on
for j = 2:freqLength-1
    if w(j)>w(j-1) && w(j)>w(j+1)
        plot(j,w(j),'r*')
        hold on
    end
end

%Remoção dos valores abaixo do limiar de um único frame

for m = 1:freqLength
    if w(m)<fitPoly(m)
        w(m) = NaN;
    end
end

%Visualização dos picos sozinhos e construção do vetor de picos

peakVector = w;
plot(w)
hold on
for j = 2:freqLength-1
    if w(j)>w(j-1) && w(j)>w(j+1)
        peakVector(j) = w(j);
        plot(j,w(j),'r*')
        hold on
    else
        peakVector(j) = NaN;
    end
end
stem(peakVector);

%teste de plot de mais de um valor por variavel horizontal
for n = 0:4
    for o=1:freqLength
        if not(isnan(peakVector(o)))
            plot(n,o,'r+')
            hold on
        end
    end
            
end

%Geração de polinômio de limiar para todos os frames
powerFramesdB = 10*log10(powerFrames);
polyOrder = 3; % ordem do polinomio
%geração de coeficientes
for i=1:timeLength
    polyCoefMatrix(i,1:polyOrder+1) = polyfit(1:freqLength,powerFramesdB(i,1:freqLength),polyOrder); 
end

%geração do polinomio
for i=1:timeLength
    fitPolyMatrix(i,1:freqLength)= polyval(polyCoefMatrix(i,1:polyOrder+1),1:freqLength); 
end

%verificação dos polinomios
for i=1:timeLength
    plot(1:freqLength,fitPolyMatrix(i,1:freqLength))
    hold on
end

%Remoção dos valores abaixo do limiar dos frames

for i=1:timeLength
    for k=1:freqLength
        if powerFramesdB(i,k)<fitPolyMatrix(i,k)
            powerFramesdB(i,k)=NaN;
        end
    end
end

%Construção da matriz de picos

peakMatrix = powerFramesdB; % inicializando para ter mesma dimensão ao fim da operação
for i=1:timeLength
    for k=2:freqLength-1
        if powerFramesdB(i,k)>powerFramesdB(i,k-1) && powerFramesdB(i,k)>powerFramesdB(i,k+1)
            peakMatrix(i,k) = powerFramesdB(i,k);
        else
            peakMatrix(i,k) = NaN;
        end
    end
end

%teste de plot de mais de um valor por variavel horizontal com ligação
%entre frequencias iguais a cada frame
for n=1:20
    for o=1:freqLength
        if not(isnan(peakMatrix(n,o)))
            plot(n,o,'bo')
            if not(isnan(peakMatrix(n,o))) && not(isnan(peakMatrix(n+1,o))) 
                line([n n+1],[o o]);
            end
            hold on
        end
    end        
end

%teste de laço que encontra a frequencia mais proxima entre n e n+1 e faz a
%ligação(sem mais nenhuma especificação como resolver duplicação)
minDistanceFreq = zeros(size(fitPoly));
freqDistances = zeros(size(fitPoly));
for n=1:100
    for o=1:freqLength-1
        if not(isnan(peakMatrix(n,o)))
            plot(n,o,'b');
            for q=1:freqLength-1
                if not(isnan(peakMatrix(n,o))) && not(isnan(peakMatrix(n+1,q)))
                    freqDistances(q)=norm([n+1 q] - [n o]);
                    if freqDistances(q+1)<freqDistances(q)
                        minDistanceFreq(q) = q+1;
                    end
                end
            end
            line([n n+1],[o minDistanceFreq(o)]);
            hold on
        end     
    end
end


% minFreqDistanceIndex = zeros;
% freqNorm = zeros(size(fitPoly));
% for n=1:10
%     for o=1:freqLength
%         if not(isnan(peakMatrix(n,o)))
%             plot(n,o,'b');
%             for q=1:freqLength-1
%                 if not(isnan(peakMatrix(n+1,q)))
%                     freqNorm(q) = norm([n+1 q]-[n o]);
%                     if freqNorm(q)<freqNorm(q+1)
%                         minFreqDistanceIndex = q;
%                     end
%                 end
%             end
%             line([n n+1],[o minFreqDistanceIndex]);
%             hold on
%         end
%     end
% end
% 
% minFreqDistanceIndex = zeros;
% freqNorm = zeros(size(fitPoly));
% for n=1:10
%     for o=1:freqLength
%         if not(isnan(peakMatrix(n,o)))
%             plot(n,o,'b');
%             for q=1:freqLength-1
%                 if not(isnan(peakMatrix(n+1,q)))
%                     freqNorm(q) = norm([n+1 q]-[n o]);
%                     if freqNorm(q)<freqNorm(q+1)
%                         minFreqDistanceIndex = q;
%                     end
%                 end
%             end
%             line([n n+1],[o minFreqDistanceIndex]);
%             hold on
%         end
%     end
% end
