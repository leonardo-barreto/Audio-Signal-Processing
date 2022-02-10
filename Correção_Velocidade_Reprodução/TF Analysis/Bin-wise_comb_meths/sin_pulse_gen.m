function [x,t] = sin_pulse_gen(f0, fs, t_i, t_f, nHarm, dacay_flag)
%SIN_PULSE_GEN Summary of this function goes here

t = 0 : 1/fs : .5 + t_f;
decay = 1:-1/((t_f-t_i)*fs):0;

x = sin(2*pi*f0*t)';
if nHarm > 1
    for harm = 2:nHarm
        factor = 1/sqrt(harm);
        x = x + factor*sin(harm*2*pi*f0*t)';
    end
end

w = ones(length(t),1);
w(1:floor(t_i*fs)) = 0;
w(floor(t_f*fs):end) = 0;

if dacay_flag
    w(floor(t_i*fs):floor(t_i*fs)+length(decay)-1) = decay;
end

x = x.*w;
    
end

