time = 0:1/48000:0.5-(1/48000);


chirp1 = 0;
for index = 1:10
    chirp1 = chirp1 + 2*chirp(time,100*index,0.5,200*index)/2/index;
end

clearvars index

spectrogram(chirp1)