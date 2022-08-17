function [spectrograms, final_freq_bins_CQT, finalTimeVec, Xcqt] = CQT_spectrogram_tensor_prep(input_sig, fs, bins_per_octave, fmax, fmin, hop)
%CQT_spectrogram_tensor_prep is a function that computes CQT spectrograms
%with different resolutions and concatenate them into a tensor
% input_sig: input signal
% fs: sampling rate
% bins_per_octave: array with the number of bins per octave
% fmax: maximum frequency of CQT
% fmin: minimum frequency of CQT
% hop: general hop size

addpath cqt_toolbox

% Computing CQT spectrograms - - - - - - - - - - - - - - - - - - - - - - 

% % Drop frequencies outside [fmin fmax]:
% w1 = 1.1*(fmin/(fs/2));
% w2 = 0.9*(fmax/(fs/2));
% [B,A] = butter(6,[w1 w2]); 
% input_sig = filtfilt(B,A,input_sig); 

% Computing rasterized complex coefficients - - - - - - - - - - - - - - 
spectrograms_all = cell(length(bins_per_octave), 3);
for ind = 1:length(bins_per_octave)
    Xcqt = cqtPerfectRast(input_sig, fmin, fmax, bins_per_octave(ind), fs, 'q', 1, 'atomHopFactor',...
                                      0.25, 'thresh', 0.0005, 'win', 'sqrt_blackmanharris');

    spec_CQT = abs(Xcqt.spCQT).^2;

    emptyHops = Xcqt.intParams.firstcenter/Xcqt.intParams.atomHOP;
    maxDrop = emptyHops*2^(Xcqt.octaveNr - 1) - emptyHops;
    droppedSamples = (maxDrop - 1)*Xcqt.intParams.atomHOP + Xcqt.intParams.firstcenter;
    TimeVec = (1:size(spec_CQT, 2))*Xcqt.intParams.atomHOP - Xcqt.intParams.preZeros + droppedSamples;

    spec_CQT(:,(TimeVec < 0) | (TimeVec > length(input_sig))) = [];
    TimeVec((TimeVec < 0) | (TimeVec > length(input_sig))) = [];

    freq_bins_CQT = Xcqt.fmin*2.^(((1:Xcqt.octaveNr*Xcqt.bins)-1)/Xcqt.bins);

    spectrograms_all{ind, 1} = spec_CQT;
    spectrograms_all{ind, 2} = TimeVec;
    spectrograms_all{ind, 3} = freq_bins_CQT;
end


% Interpolating CQTs - - - - - - - - - - - - - - - - - - - - - - - - - - 
finalTimeVec = 0:hop:TimeVec(end)/fs;

final_freq_bins_CQT = spectrograms_all{end, 3};

% N_w_str = '';
for CQT_ind = 1:length(bins_per_octave)
%     N_w_str = [N_w_str sprintf('_%d', bins_per_octave(CQT_ind))];
    [X,Y] = meshgrid(spectrograms_all{CQT_ind, 2}/fs, spectrograms_all{CQT_ind, 3});
    [Xq, Yq] = meshgrid(finalTimeVec, final_freq_bins_CQT);

    spectrograms_all{CQT_ind, 1} = interp2(X, Y, spectrograms_all{CQT_ind, 1}, Xq, Yq, 'linear', 0);

%     if plot_flags(1)
%         plot_spec_CQT(sprintf('CQT_%dbins', bins_per_octave(CQT_ind)), compress_dB_norm(spectrograms_all{CQT_ind, 1},80), ...
%                                            finalTimeVec, final_freq_bins_CQT, Xcqt, plot_par, figs_path);
%     end
end

spectrograms = zeros(length(final_freq_bins_CQT), length(finalTimeVec), length(bins_per_octave));
for ind = 1:length(bins_per_octave)
    spectrograms(:,:,ind) = spectrograms_all{ind, 1};
end


% Normalizing CQTs - - - - - - - - - - - - - - - - - - - - - - - - - - -  
orig_energy = sum(input_sig.^2);

for ind = 1:size(spectrograms,3)
    spec_energy = sum(sum(spectrograms(:,:,ind)));
    spectrograms(:,:,ind) = spectrograms(:,:,ind)*orig_energy/spec_energy;
end

end