function alphasEST_result = filter_alphas(alphasEST, ord)

% Separating the alphas
% alphasEST_diff = diff(alphasEST);

% figure;
% plot(alphasEST(:,1),'r'); hold on;
% plot(alphasEST(:,2),'k');
% plot(alphasEST_diff(:,1), 'b');
% plot(alphasEST_diff(:,2), 'g');
% legend('Est. Alphas 1', 'Est. Alphas 2', 'Diff. Est. Alphas 1', 'Diff. Est. Alphas 2')

% lpfilter = ones(1,ord)/ord;
lpfilter = hamming(ord)/sum(hamming(ord));

alphasEST_result = zeros(size(alphasEST));

for ind = 1:size(alphasEST,2)
    alphasEST_result(:,ind) = filtfilt(lpfilter, 1, alphasEST(:,ind));
end
end