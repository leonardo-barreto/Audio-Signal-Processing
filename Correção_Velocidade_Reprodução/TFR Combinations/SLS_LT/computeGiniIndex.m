function sparsity = computeGiniIndex(S)

sparsity = zeros(1,size(S,3));

for tfrInd = 1:size(S,3)
    if sum(sum(S(:,:,tfrInd))) > 0
        vector = reshape(S(:,:,tfrInd), 1, size(S(:,:,tfrInd),1)*size(S(:,:,tfrInd),2));
        sortedVector = sort(vector);
        sortedVector = sortedVector/sum(sortedVector);

        N = length(sortedVector);
        k = 1:N;
        auxVector = (N - k + 0.5)./N;

        sparsity(tfrInd) = 1 - 2*sum(sortedVector.*auxVector);
    else
        sparsity(tfrInd) = 0;
    end
end

end