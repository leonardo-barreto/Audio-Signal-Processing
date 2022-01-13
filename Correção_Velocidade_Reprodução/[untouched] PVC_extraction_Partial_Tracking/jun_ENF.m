function [track_F_mean]=jun_ENF(Partials, time, fs, num_tracks)
% Track ENF in partials
% 
%
t_min=0.5; % minimun track length
t_min_frames=t_min*fs/(time(2)-time(1)); % minimun number of frames for track length
M = length(time);
track_F = -inf+zeros(M,num_tracks);
track_A = -inf+zeros(M,num_tracks);


for m = 1:M    
    active = Partials{m}(:,4);
    track_F(m,active) = Partials{m}(:,1);
    track_A(m,active) = Partials{m}(:,2);
end

track_index=sum(track_F~=-inf,1)>t_min_frames;
track_F_filtered=track_F(:,track_index);
track_F_filtered(track_F_filtered==-Inf)=NaN; %convert -Inf to NaN to calculate mean
track_F_mean=mean(track_F_filtered,'omitnan');
track_F_filtered=track_F(:,track_index);
track_F_normalized=track_F_filtered./track_F_mean;

track_F_normalized(track_F_normalized==-Inf)=NaN; %convert -Inf to NaN to mean

%track_F_mean=mean(track_F_normalized,2,'omitnan');
variance_tracks_clustering=movvar(track_F_normalized,10)>2e-5;
clustering=sum(movvar(track_F_normalized,30)>2e-5);
track_F_mean=mean(track_F_normalized(:,clustering==0),2,'omitnan');

figure
plot(time/fs, fs/(2*pi)*track_F_normalized(:,clustering==0), 'linewidth',1.5);

end