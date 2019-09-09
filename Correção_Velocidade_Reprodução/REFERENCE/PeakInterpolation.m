function [iploc, ipmag] = PeakInterpolation (mX, ploc)
mXdb = 20*log10(mX);        % converting the magnitude spectrum into dB
val = mXdb(ploc);             % peak amplitudes in dB
lval = mXdb(ploc - 1);          % amplitudes on peaks' left
rval = mXdb(ploc + 1);          % amplitudes on peaks' right
iploc = ploc + 0.5*(lval - rval)./(lval - 2*val + rval);    % estimated peaks' locations
ipmag = db2mag(val - 0.25*(lval - rval));           % estimated peaks' magnitude in dB
end