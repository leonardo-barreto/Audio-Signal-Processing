function alphasEST_result = filter_alphas_median(alphasEST, ord)

alphasEST_result = alphasEST;

% Only odd orders will be considered, even orders will be reduced
if ~mod(ord,2); ord = ord - 1; end;

for sourceInd = 1:size(alphasEST,2)
    for sampleInd = (ord+1)/2 : size(alphasEST,1) - (ord-1)/2
        vector = alphasEST(sampleInd - (ord-1)/2: sampleInd + (ord-1)/2, sourceInd);
        vectorOrd = sort(vector);
        alphasEST_result(sampleInd, sourceInd) = vectorOrd((ord+1)/2);
    end
end

end