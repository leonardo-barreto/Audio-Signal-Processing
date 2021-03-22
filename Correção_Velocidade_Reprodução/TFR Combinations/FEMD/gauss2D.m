function G = gauss2D(N,sigma)

x  = -(N-1)/2:(N-1)/2;
gb = 1/sqrt(2*pi*sigma(1)^2)*exp(-x.^2/(2*sigma(1)^2)); gb = gb(:);
gk = 1/sqrt(2*pi*sigma(2)^2)*exp(-x.^2/(2*sigma(2)^2)); gk = gk(:);

G = gb*(gk');

end

