
% Generating the 2D windows
size_w = [33 33*2];
wr = window(@hamming, size_w(1));
wc = window(@hamming, size_w(2));

wc(ceil(end/2) + 1 : end) = 0;
[maskr, maskc] = meshgrid(wc,wr);
W = maskc.*maskr;
W(:,size_w(2)/2:size_w(2)/2) = [];

offset_x = (size_w(2)-1)/2;
offset_y = (size_w(1)-1)/2;

x_ticks = -(size_w(2)-1)/2 + offset_x : 0 + offset_x;
y_ticks = -(size_w(1)-1)/2 + offset_y : (size_w(1)-1)/2 + offset_y;

figure;
black = abs(1-gray);
imagesc(x_ticks, y_ticks, W)
colormap(black)
colorbar
% xlim([-8,1])
xlabel('m')
ylabel('k')
set(gca,'YDir','normal');

return

% Plotting samples from long window
windowSizeFactor = 4;
for freq = 0:1/(4*windowSizeFactor):1-1/(4*windowSizeFactor)
    hold on;
    plot(1:windowSizeFactor:32, freq*ones(32/windowSizeFactor), 'ro');
end

% Plotting samples from medium window
windowSizeFactor = 2;
for freq = 0:1/(4*windowSizeFactor):1-1/(4*windowSizeFactor)
    plot(1:windowSizeFactor:32, freq*ones(32/windowSizeFactor), 'bo', 'MarkerSize', 10);
end

% Plotting samples from short window
windowSizeFactor = 1;
for freq = 0:1/(4*windowSizeFactor):1-1/(4*windowSizeFactor)
    plot(1:windowSizeFactor:32, freq*ones(32/windowSizeFactor), 'ko', 'MarkerSize', 5);
end
hold off;


