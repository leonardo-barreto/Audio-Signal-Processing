function [ x_out ] = make_mono( x_in )
%MAKE_MONO Summary of this function goes here
%   Detailed explanation goes here
if size(x_in,2) > 1
    x_out = mean(x_in,2);
end

end

