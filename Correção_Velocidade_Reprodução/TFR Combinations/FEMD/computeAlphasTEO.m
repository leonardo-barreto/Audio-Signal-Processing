function [alphasTEO] = computeAlphasTEO(x, fs, ffund, N, hop)
% Compute theoretical values of alpha for synthetic signals

noverlap   = N - hop; 
Num_frames = fix((length(x)-noverlap)/(N-noverlap)); 
alphasTEO   = zeros(Num_frames,1);
for i=1:Num_frames
    % Dividing signal in frames:
    t0 = (i-1)*hop+1;
    fn = ffund(t0:t0+N-1);
    % Estimating alpha value:
    auxn  = 0:(N-1); T = (N-1)/2/fs;
    tempf = auxn/fs - T;
    % Finding central frequency:
    fcentral = interp1(tempf,fn(1,:),0);
    frame    = fn(1,:)/fcentral;
    % Compute and keep coefficients:
    coefs = polyfit(tempf,frame,1);
    parametro = coefs(1)/coefs(2);
    alphasTEO(i) = parametro;
end

end

