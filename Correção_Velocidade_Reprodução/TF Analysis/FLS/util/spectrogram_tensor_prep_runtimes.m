function [spectrograms, Freq, Time, runtimes ] = spectrogram_tensor_prep_runtimes(input_sig, fs, N_w, NFFT, overlap_short)
%SPECTROGRAM_COMB is a function that combines different resolution
%spectrograms by a point-wise weighted geometric average.
% input_sig: input signal
% fs: sampling rate
% sig_name: string which will be used to name the figures generated
% N_w: vector containing the window lentgh used for each spectrogram computation
% Beta: vector containing the chosen beta parameters
% plot_par: vector containing the time em frequency limits of the plots and
% the font size [t_min, t_max, f_min, f_max, font_size]
   %            ex: plot_par = [1, 5, 0, 5000,25]; 
% plot_flags: flags used to choose which spectrograms will be plotted 
%               [spec1, spec2, spec3, combined]
% redLines: vector containing the instants in which will be plotted red
% dashed vertical lines

% Computing spectrograms - - - -
spectrogramsBeforeTrim = cell(length(N_w));

runtimes = zeros(size(N_w));

N_w = sort(N_w); % just to make shure they are ordered

tic
[spec_short, Freq, Time] = spectrogram(input_sig, hanning(N_w(1)), overlap_short, NFFT, fs);
spectrogramsBeforeTrim{1} = spec_short;
runtimes(1) = toc;

for ind = 2:length(N_w)
    x_shifted = [zeros((N_w(ind) - N_w(1))/2,1); input_sig];
    tic
    spectrogramsBeforeTrim{ind} = spectrogram(x_shifted, hanning(N_w(ind)), N_w(ind) - N_w(1) + overlap_short, NFFT, fs);
    runtimes(ind) = toc;
end

% Triming the lengths - - - - - - - -

sizes = zeros(1,length(N_w));

for ind = 1:length(N_w)
    sizes(ind) = size(spectrogramsBeforeTrim{ind} ,2);
end

spec_length = min(sizes);

spectrograms = zeros(size(spec_short,1), spec_length, length(N_w));

orig_energy = sum(input_sig.^2);

for ind = 1:length(N_w)
    spectrograms(:,:,ind) = abs(spectrogramsBeforeTrim{ind}(:, 1:spec_length)).^2;
    spec_energy = sum(sum(spectrograms(:,:,ind)));
    spectrograms(:,:,ind) = spectrograms(:,:,ind)*orig_energy/spec_energy;
end

Time = Time(1:spec_length);

end