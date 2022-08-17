function plot_window_S(size_w_S, name, font_size)
% PLOT windows
size_w_S = ceil(size_w_S);

if mod(size_w_S(1), 2) == 0    
    size_w_S(1) = size_w_S(1) + 1;
    fprintf('WARNING: size_w_S(1) must be an odd number! Using %d instead!\n', size_w_S(1));
end

if mod(size_w_S(2),2) == 0
    size_w_S(2) = size_w_S(2) + 1;
    fprintf('WARNING: size_w_S(2) must be an odd number! Using %d instead!\n', size_w_S(2));
end

% Generating the 2D windows
wr = window(@hamming, size_w_S(1));
wc = window(@hamming, size_w_S(2));
[maskr, maskc] = meshgrid(wc,wr);
W = maskc.*maskr;
x_ticks = -(size_w_S(2)-1)/2 : (size_w_S(2)-1)/2;
y_ticks = -(size_w_S(1)-1)/2 : (size_w_S(1)-1)/2;

figure;
black = abs(1-gray);
imagesc(x_ticks, y_ticks, W)
colormap(black)
% colorbar
xlabel('m')
ylabel('k')

figProp = struct('size', font_size, 'font','Times','lineWidth',2,'figDim',[1 1 500 300]);
figFileName = name;
formatFig(gcf,figFileName,'en',figProp);

end