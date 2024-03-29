function w = hanningz(n)
%HANNING HANNINGZ(N) returns the N-point Hanning window in a column vector.

%	Copyright (c) 1984-94 by The MathWorks, Inc.
%       $Revision: 1.4 $  $Date: 1994/01/25 17:59:15 $

w = .5*(1 - cos(2*pi*(0:n-1)'/(n)));