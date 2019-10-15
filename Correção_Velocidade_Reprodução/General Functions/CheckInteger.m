%Funcao para checar se a entrada é um inteiro

%Saidas: True/False

    function output = CheckInteger(input)
        output = (mod(input,1) == 0);
    end