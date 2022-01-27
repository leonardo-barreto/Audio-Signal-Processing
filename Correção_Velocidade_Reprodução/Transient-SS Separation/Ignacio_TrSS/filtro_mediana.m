function [ ModuloTo, ModuloIm ] = filtro_mediana(modulo,num_tonal,num_percusivo)
%filtro_mediana 
%% Input
% modulo Spectrogram module matrix
% num_tonal filter length of steady-state median filter 
% num_percusivo filter length of transient median filter
%% Output
% ModuloTo Steady-state components spectrogram
% Modulo_Im Transient components spectrogram


ModuloTo=zeros(size(modulo));
ModuloIm=zeros(size(modulo));
%ModH=zeros(size(Modulo));
%ModRes=zeros(size(Modulo));
 %Mediana en tiempo S_h
%h=waitbar(0,'Filtro de Mediana armónico');
largo=size(modulo,2);
textprogressbar('Filtering: ')
    
for  k=1:largo
    
    
        if (k>num_tonal)&&(k<=size(modulo,2)-num_tonal)
            ModuloTo(:,k)=median(modulo(:,k-num_tonal:k+num_tonal),2);
        else
            if k<=num_tonal
            ModuloTo(:,k)=median(modulo(:,1:k+num_tonal),2);
            else
            ModuloTo(:,k)=median(modulo(:,k-num_tonal:k),2);    
            end
            
        end
        textprogressbar((k/2/largo)*100);
        % waitbar(k/largo);
end
%close(h);
% Mediana en Frecuencia, S_f
ancho=size(modulo,1);
%h=waitbar(0,'Filtro de Mediana impulsivo');
    
    for k=1:ancho
        if (k>num_percusivo)&&(k<=size(modulo,1)-num_percusivo)
            ModuloIm(k,:)=median(modulo(k-num_percusivo:k+num_percusivo,:),1);
        else
            if k<=num_percusivo
            ModuloIm(k,:)=median(modulo(1:k+num_percusivo,:),1);
            else
                
                ModuloIm(k,:)=median(modulo(k-num_percusivo:k,:),1);
            end
            %ModuloIm(k,:)=0;
        end
 %   waitbar(k/ancho);
    textprogressbar((1/2+k/2/ancho)*100);
    end
    
%close(h);


ToMask=ModuloTo.^2./(ModuloTo.^2+ModuloIm.^2);
TrMask=ModuloIm.^2./(ModuloTo.^2+ModuloIm.^2);
ModuloTo=modulo.*ToMask;
ModuloIm=modulo.*TrMask;
textprogressbar('Done')
end

