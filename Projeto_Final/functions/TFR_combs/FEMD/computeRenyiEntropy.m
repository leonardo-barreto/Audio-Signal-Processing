function renyiEntropy = computeRenyiEntropy(S, Nf, hop, renyiOrd)

Delta = log2(0.5*hop/Nf);

renyiEntropy = zeros(1,size(S,3));

for tfrInd = 1:size(S,3)
    globalSum = sum(sum(S(:,:,tfrInd)));
    renyiEntropy(tfrInd) = (1/(1-renyiOrd))*log2((sum(sum((S(:,:,tfrInd)/globalSum).^renyiOrd)))) + Delta;
end

end