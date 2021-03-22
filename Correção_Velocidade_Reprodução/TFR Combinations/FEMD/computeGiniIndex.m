function sparsity = computeGiniIndex(S)

sparsity = zeros(1,size(S,3));

for TFR_Ind = 1:size(S,3)
    if sum(sum(S(:,:,TFR_Ind))) > 0
        vector = reshape(S(:,:,TFR_Ind), 1, size(S(:,:,TFR_Ind),1)*size(S(:,:,TFR_Ind),2));
        sortedVector = sort(vector);
        sortedVector = sortedVector/sum(sortedVector);

        N = length(sortedVector);
        k = 1:N;
        auxVector = (N - k + 0.5)./N;

        sparsity(TFR_Ind) = 1 - 2*sum(sortedVector.*auxVector);    
    else
        sparsity(TFR_Ind) = 0;
    end
end

end