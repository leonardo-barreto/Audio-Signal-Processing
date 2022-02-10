function alphasEst_result = find_alphas(xpdf, ypdf, percent_trsh, NumOfSources)

[peaks, locs] = findpeaks(ypdf);
locs = xpdf(locs);

[peaks_sorted, inds_sorted]= sort(peaks);

if ~isempty(peaks_sorted)
    inds_sorted(peaks_sorted < percent_trsh*peaks_sorted(end)) = inds_sorted(end);
    inds = inds_sorted;

    locs_sorted = locs(inds);
    locs_sorted = locs_sorted(end:-1:1);

    if length(locs_sorted) < NumOfSources
        alphasEst_result = [locs_sorted ; zeros(NumOfSources - length(locs_sorted),1)];
    else
        alphasEst_result = locs_sorted(1:NumOfSources)';
    end

    alphasEst_result = sort(alphasEst_result);
else
    alphasEst_result = zeros(1,NumOfSources);
end

end