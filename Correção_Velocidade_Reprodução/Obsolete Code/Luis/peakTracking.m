function [tfreqn, memoryOut] = peakTracking (pmag, pfreq, tfreq, freqDevOffset, freqDevSlope, memoryIn)

tfreqn = zeros(1, length(tfreq));           % inicializa array com as frequencias da saída
pindexes = find(pfreq);            % índices dos picos atuais
incomingTracks = find(tfreq);               % índices da trilha de entrada
newTracks = zeros(1, length(tfreq))-1;      % inicializa novas trilhas com -1
[~, magOrder] = sort(pmag(pindexes),'descend');
pfreqt = pfreq;                    % copia os picos atuais para um array temporario
% memoryOut = memoryIn;
%% Continuação das trilhas já existentes
if ~isempty(incomingTracks)                 % se existir trilhas de entrada
    for i = 1:length(magOrder)
        if isempty(incomingTracks)          % sai do loop quando todas as trilhas já existentes acharem uma continuacao
            break
        end
        [~, track] = min(abs(pfreqt(magOrder(i)) - tfreq(incomingTracks)));          % frequencia mais proxima do pico
        freqDistance = abs(pfreqt(magOrder(i)) - tfreq(incomingTracks(track))); % criterio da distancia minima
        if freqDistance < (freqDevOffset + freqDevSlope*pfreqt(magOrder(i)))        % tolerancia
            newTracks(incomingTracks(track)) = magOrder(i);                     % atribui o indice do pico ao indice da trilhas
            incomingTracks(track) = [];                                         % remove a trilha atribuida das trilhas de entrada
        end
    end
end
indext = find(newTracks ~= -1);             % indices das trilhas assimiladas com frequencias do frame atual
if ~isempty(indext)                         % assimilando as trilhas continuadas ao vetor tfreqn
    indexp = newTracks(indext);             % indices dos picos assimilados(newTracks é um vetor com os indices dos picos)
    tfreqn(indext) = pfreqt(indexp);        % atribui frequencias de saida provenientes de continuacao de trilhas (e nao a partir de novas trilhas)
    pfreqt(indexp) = [];                    % depois de atribuir, remove os picos assimilados
end

memoryIn(intersect(find(memoryIn),find(tfreqn))) = 0;    % zera os picos virtuais que encontraram um pico
memoryIn(intersect(find(tfreq),find(~tfreqn))) = memoryIn(intersect(find(tfreq),find(~tfreqn))) + 1;  % incrementa o contador ao copiar o pico 
tfreqn(intersect(find(tfreq),find(~tfreqn))) = tfreq(intersect(find(tfreq),find(~tfreqn))); % copia os picos do frame antigo que nao encontraram candidatos para o frame atual
memoryOut = memoryIn;   % copia o vetor de memoria atualizado para a memoria de saida

%% Criar novas trilhas à partir de picos não utilizados
emptyt = find(~tfreq);                      % indices das trilhas de entrada vazias
% peaksleft = find(pfreqt);                   % indices dos picos ainda nao assimilados que sobraram (variante para picos ordenados por magnitude)
if (~isempty(pfreqt) && (length(emptyt) >= length(pfreqt)))      % caso existam MAIS espacos vazios do que picos nao utilizados
    tfreqn(emptyt(1:length(pfreqt))) = pfreqt;                          % completa as trilhas vazias com os picos nao assimilados
    %tfreqn(emptyt(1:length(peaksleft))) = pfreqt(peaksleft);            % variante da linha de código acima
elseif (~isempty(pfreqt) && (length(emptyt) < length(pfreqt)))   % caso existam MENOS espacos vazios do que picos nao utilizados
    tfreqn(emptyt) = pfreqt(1:length(emptyt));
    tfreqn = [tfreqn pfreqt(1+length(emptyt):end)];                     % cria novas trilhas com os picos que sobraram
    %tfreqn(emptyt) = pfreqt(peaksleft(1:length(emptyt)));               % variante da linha de código acima
    %tfreqn = [tfreqn pfreqt(peaksleft(1+length(emptyt):end))];          % variante da linha de código acima
end

end