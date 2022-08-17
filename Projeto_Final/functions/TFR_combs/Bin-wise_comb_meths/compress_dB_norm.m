function X_comp = compress_dB_norm(X, range)

X_comp = X/max(max(max(X)));
X_comp = 20*log10(X_comp)*(1/range) + 1;
X_comp(X_comp<0) = 0;

end