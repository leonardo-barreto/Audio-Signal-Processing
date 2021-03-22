function [STFChT, freq, time] = fanchirpogram(varargin)

x = varargin{1};
w = varargin{2};
hop = varargin{3};
fs = varargin{4};
alphas = varargin{5};

if length(varargin) < 6
    N_FFT = length(w);
else
    N_FFT = varargin{6};
end

K = length(w);

time = K/(2*fs) : hop/fs : ((floor((length(x) - K) / hop) + 1)*hop + K/2 - 1)/fs;
freq = (0:K/2)*fs/K;

STFChT = zeros(K/2+1, length(time));

% warping functions
psi_warp = @(x, alfa) (-1+sqrt(1+2*alfa*x))/alfa;
phi_warp = @(x, alfa) (1+1/2*alfa*x)*x;

if length(alphas) < length(time)
    alphas = alphas(1);
    
    M = 2^nextpow2(ceil(K / (1-abs(alphas/fs) * K/2)));
    if M > K
        w = resample(w, M, K);
    end
    
    % Resampling:
    if alphas ~= 0 
        tn = ((0:1:K-1)-(K-1)/2)/fs;
        tr = phi_warp(tn(1),alphas) + ...
                            ((0:1:M-1)+1/2)*(tn(end)-tn(1))/M;
        tt = psi_warp(tr,alphas);
    else
        tn = ((0:1:K-1)-(K-1)/2)/fs;
        tt = tn(1) + ((0:1:M-1)+1/2)*(tn(end)-tn(1))/M;
    end

    frames = zeros(length(w), length(time));
    for frame_ind = 1:length(time)
        t0 = (frame_ind-1)*hop + 1;
        xn = x(t0 : t0+K - 1);
        xr = interp1(tn,xn,tt);
        frames(:, frame_ind) = xr.*w';
    end
    
    X = fft(frames, N_FFT*M/K, 1);
    
    STFChT = X(1 : N_FFT/2+1, :);
else
    for frame_ind = 1:length(time)
        alpha_max = max(alphas);
        M = 2^nextpow2(ceil(K / (1-abs(alpha_max/fs) * K/2)));
        w = resample(w, M, N);

        % Resampling:
        if alphas(frame_ind) ~= 0 
            tn = ((0:1:K-1)-(K-1)/2)/fs;
            tr = phi_warp(tn(1),alphas(frame_ind)) + ...
                                ((0:1:M-1)+1/2)*(tn(end)-tn(1))/M;
            tt = psi_warp(tr,alphas(frame_ind));
        else
            tn = ((0:1:K-1)-(K-1)/2)/fs;
            tt = tn(1) + ((0:1:M-1)+1/2)*(tn(end)-tn(1))/M;
        end
        
        t0 = (frame_ind-1)*hop + 1;
        xn = x(t0 : t0+K - 1);
        xr = interp1(tn,xn,tt);
        Xr = fft(xr'.*w, N_FFT*M/K);
        STFChT(:,frame_ind) = Xr(1 : N_FFT/2+1);
    end
end
end
