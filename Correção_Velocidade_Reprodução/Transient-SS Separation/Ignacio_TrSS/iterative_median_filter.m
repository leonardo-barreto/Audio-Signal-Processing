function [ ModuloTo,ModuloRes, ModuloImp ] = iterative_median_filter(modulo,num_ss,num_t,niter)
%%function [ ModuloTo,ModuloRes, ModuloImp ] = filtro_mediana_tricomponente(modulo,num_ss,num_t,niter)
% Iterative filtering of the spectrogram

%%Initializing variables
ModuloTo=zeros(size(modulo));
ModuloImp=zeros(size(modulo));
ModuloRes=zeros(size(modulo));
ModuloToTemp=zeros(size(modulo));
ModuloImpTemp=zeros(size(modulo));

[ModuloTo,ModuloImp]=filtro_mediana(modulo,num_ss,num_t);%First processing
%Iterative procesing
if (niter~=0)
    for i=1:niter
        [ModuloTo,ModuloImpTemp]=filtro_mediana(ModuloTo,num_ss,num_t);
        ModuloRes=ModuloRes+ModuloImpTemp;
        [ModuloToTemp,ModuloImp]=filtro_mediana(ModuloImp,num_ss,num_t);
        ModuloRes=ModuloRes+ModuloToTemp;
    end
end 

end


