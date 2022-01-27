function phase_out=princarg(phase_in)
% princarg, function to compute the principal phase argument in the
% range [-pi, pi].
%
% usage: phase_out = princarg(phase_in)

phase_out=mod(phase_in+pi,-2*pi)+pi;
