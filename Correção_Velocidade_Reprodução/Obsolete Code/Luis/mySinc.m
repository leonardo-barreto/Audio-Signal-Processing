function  output  = mySinc(x)
% Funcao que retorna o valor de sen(x)/x. X em radianos.

output = ones(size(x));
nonzero = x~=0;
output(nonzero) = sin(x(nonzero))./x(nonzero);

% if x == 0
%     output = 1;
% else
%     output = sin(x)./x;
% end
end

