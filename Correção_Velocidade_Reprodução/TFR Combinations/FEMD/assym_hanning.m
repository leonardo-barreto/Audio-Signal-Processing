function w = assym_hanning(N1,N2)
w1 = hamming(N1); 
w2 = hamming(N2);
w = [w1(1:end/2) ; w2(end/2 + 1:end) ; zeros(N1/2 - N2/2 , 1)];
end