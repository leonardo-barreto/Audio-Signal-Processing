function [xtfreq trackAverage timeAxis freqAxis blockAxis frame samples] = SinusoidalAnalysis (xSound, fs, N, M, freqDevSlope, minTrackLength, memoryLength, offsetTPSW)
%% CONSTANTES
R = (M-1)/2;                                % fator de overlap
maxnSines = 60;                            % numero maximo de senoides em um frame
freqDevOffset = 0;                         % desvio minimo em 0Hz

%% ANÁLISE SENOIDAL
frame = fix((length(xSound) - M)/R);                % numero de frames
% disp(frame);
samples = fix((0:frame-1)*R + (M+1)/2);              % numero da amostra instantanea referente a cada quadro
timeAxis = samples/fs;                              % instante de tempo referente a cada quadro
freqAxis = (0:N/2)*fs/N;                            % escalando o eixo x para exibir as frequencias
blockAxis = 1:frame;
tfreq = [];
memory = zeros(1,maxnSines);
for i = 1:frame
    x1 = xSound(((i-1)*R)+1:((i-1)*R + M));         % seleciona o frame
    % TRANSFORMADA DE FOURIER
    [mX mXdb] = dft(x1,N);
    % DETECÇÃO DE PICOS
    ploc = peakDetect(mX,offsetTPSW);
    % INTERPOLAÇÃO
    [iploc ipmag] = peakInterpolation (mXdb, ploc);
    ipfreq = (iploc)*fs/N;
    % RASTREAMENTO DE PICOS
    [tfreq, memory] = peakTracking (ipmag, ipfreq, tfreq, freqDevOffset, freqDevSlope, memory);
    tfreq = tfreq(1:min(maxnSines, length(tfreq)));              % limitar o numero maximo de trilhas em um frame para maxnSines
    jtfreq = zeros(1,maxnSines);                                    % vetor temporário de saída
    jtfreq(1:length(tfreq)) = tfreq;                              % guardando as frequencias das trilhas no vetor temporário
    if i == 1                               % caso seja o primeiro frame, inicializa a saída
        xtfreq = jtfreq;
    else                                    % para o resto, adiciona as frequencias do frame atual na saída
        xtfreq = [xtfreq ; jtfreq];
    end
    if i>memoryLength
        memdel = memory>memoryLength;             % indices das trilhas que estouraram a memoria
        xtfreq(end-memoryLength:end,memdel) = 0;        % deleta os picos virtuais criados
        memory(memory>memoryLength) = 0;                % zera o contador das trilhas que estouraram a memoria
    end
%     disp(floor(100*i/frame));clc;
end
% REMOÇÃO DE TRILHAS CURTAS
xtfreq = removingShortTracks(xtfreq, minTrackLength);
% PLOTS DOS RESULTADOS
trackAverage = averageTracks(xtfreq);

% CORREÇÃO DOS DESVIOS
% trackAverage(isnan(trackAverage)) = 1;
% resampleOutput = TimeVaryingResample(xSound', trackAverage, samples);
% wavwrite(resampleOutput, fs, 'paulistanaTrechoCorrigida');

xtfreq(xtfreq==0) = NaN;
% figure; plot(xtfreq);
end