function [Ptot,Npicos] = estimaPicos(X, mX, fs, usaSSE)
% Authors: J. Bonada, X. Serra, X. Amatriain, A. Loscos
% Estima os picos 
% X:  matriz de transformadas para cada frame
% mX: módulo da transformada para todas as frequências (para cada frame)
% fs: frequência de amostragem do sinal de entrada
% usaSST: se =1, utiliza SSE para remover o piso de ruído do espectro
% Ptot: picos estimados para cada frame
% Npicos: número de picos de cada frame
%--------------------------------------------------------------------------
Nframes   = size(mX,1);
Namostras = size(mX,2);
N2   = Namostras/2+1;
Nfft = Namostras;
Ptot = cell(1,Nframes); 
% Estimando os picos para cada frame:
for m=1:Nframes
  % Spectrum (only positive frequencies):
  spec = X(m,1:(Nfft/2+1)); %plot(abs(spec)); pause
  % Unwrapped phase spectrum
  pX = unwrap(angle(spec));             
  % abs(spec) for all frequencies:
  S  = mX(m,:); 
  % Espectro em dB:
  mXorig = 20*log10(S(:));
  % Retira o piso de ruído?
  if usaSSE==1
      % Aqui vai o pré-processamento do artigo [2]:------------------------
      % This method should be used with the whole spectrum!
      E = calculaSSE(S); %estima o piso de ruído
      E = 20*log10(E(:));
      % Retirando o piso de ruído:  
      mXaux  = mXorig(:)-E(:); % espectro sem o piso de ruído (em dB)
      %--------------------------------------------------------------------
  else
      mXaux = mXorig(:);
  end
  % Trabalhando só com as frequências positivas do espectro:
  mXorig = mXorig(1:N2); % espectro original (em dB) 
  mXaux  = mXaux(1:N2);  % espectro sem o piso de ruído (em dB)    
  % Determinando o limiar:
  t2 = min(mXaux); 
  t1 = max(mXaux)-t2;
  %t0 = max(mXaux); t = t0;
  nit = 0; t=t1+t2; 
  npicos = 0; %contador = 0;
  while (npicos<=80)%&&(contador<=10)
      % Acha os picos:
      ploc = 1 + find( (mXaux(2:(N2-1))>t).*(mXaux(2:(N2-1))>mXaux(3:N2)) ...
               .*(mXaux(2:(N2-1))>mXaux(1:(N2-2))) );
      % Atualiza o número de picos:       
      npicos = length(ploc);
      % Incrementa o limiar:
      alpha = 1/100*abs(log(2*t1/(t1-t2)));
      nit   = nit+1;
      t     = t1*exp(-alpha*nit)+t2;
      % Incrementa o contador:
      %contador = contador+1;
  end  
  clc
  % ERRO NO CÓDIGO AQUI?===================================================
  % São picos no espectro original?
  delPicos = [];
  for pico=1:npicos
      bin_atual = ploc(pico);   p2 = mXorig(bin_atual);
      bin_ant   = ploc(pico)-1; p1 = mXorig(bin_ant);
      bin_post  = ploc(pico)+1; p3 = mXorig(bin_post);
      if (p2>p1)&&(p2>p3)
      else
          delPicos = [delPicos pico];
      end
  end
  % Deleta os picos encontrados que não são picos:-------------------------
  ploc(delPicos) = [];
  npicos = length(ploc);
  %========================================================================
  % Limitando o numero maximo de picos para 80 (retirar os de menor 
  % amplitude após "planificar" o espectro!):  
  P = 80;
  if npicos>P
      [~,I] = sort(mXaux(ploc));
      I = I(end:-1:1);
      I = sort(I(1:P));
      ploc = ploc(I);
  end
  % Refine peak values: 
  [iploc,ipmag] = peakinterp(mXorig,pX,ploc);
  % iploc e ipmag tem que ter tamanho=P:
  if npicos<P
      iploc = [iploc; ((npicos+1):P).'];
      ipmag = [ipmag; mXorig(((npicos+1):P))]; 
  end
  % Mapeando as posições dos bins em frequência (Hz):
  %pfreq = (iploc-1)*fs/Nfft;
  pfreq = iploc;
  % Checando resultados:---------------------------------------------------
  % Escala de frequências da DFT:
%   fdft  = ((1:N2)-1)*fs/Nfft;
%   if m>=80
%       figure(1);
%       %plot(fdft,mXorig,'linewidth',1.5)
%       plot(fdft,mXaux,'linewidth',1.5)
%       hold on
%       %plot(pfreq,ipmag,'rx','linewidth',1.5);
%       plot((ploc-1)*fs/Nfft,mXaux(ploc),'rx','linewidth',1.5);
%       plot(fdft,t*ones(size(fdft)),'--k','linewidth',1.5)
%       hold off
%       xlabel('Frequency (Hz)','FontName','Times','FontSize',16);
%       %ylabel('Amplitude (dB)','FontName','Times','FontSize',16)
%       title('Peak Detection w/ Modification I','FontName','Times','FontSize',16)
%       xlim([0 fdft(end)]);% ylim([-150 0])
%       grid
%       pause
%   end
  %------------------------------------------------------------------------
  % Criando o vetor Pm:
  Npicos = length(iploc);
  Pm(1:Npicos,1) = pfreq;
  Pm(1:Npicos,2) = db2mag(ipmag);
  % Adicionando a Ptot:
  Ptot(1,m) = {Pm};
end

