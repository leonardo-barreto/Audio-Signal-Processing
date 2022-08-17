function [d, t, f, a] = harmonicsource(fs, dur, nharm, f0, varargin)
% Roughly simulates an harmonic sound source
% input:
% 	fs    - sampling frequency
% 	dur   - duration in seconds
%	f0    - fundamental frequency: 
%		f0(1) = initial f0
%		f0(2) = final  f0 (can be omitted if no evolution in f0 is desired)
%	nharm - number of harmonics
%
% additional options to be controlled (defaults in parenthesis): 
%       'verbose' (0)             : 0 - no messages, 1 - warnings and informative messages, 2 - plots       
%        'f0type' ('constant')    : type of function that models the evolution of f0 in time, from f0(1) up to f0(2)
%                                   other options : 'linear', 'quadratic', 'sinusoidal'
%    'amplitudes' (ones(nharm,1)) : vector of relative amplitudes of harmonics (they are normalized by default)
% 'inharmonicity' (zero)          : inharmonicity factor, used as fn = n * f0  * sqrt(1 + (n^2-1) * inharmonicity)
%                                   based on A. P. Klapuri.  Multiple fundamental frequency estimation based on 
%                                                            harmonicity and spectral smoothness. IEEE Trans. Speech 
%                                                            and Audio Proc., 11(6), 804?816, 2003
%       'amptype' ('constant')    : type of function that models the amplitude evolution of the harmonics
%                                   other options: 
%                                                 'attack-only' - amplitude grows monotonously 
%                                                  'decay-only' - amplitude falls monotonously 
%                                                'attack-decay' - amplitude grows monotonously during the attack time 
%                                                                 then falls monotonously the rest of the time
%                                                        'adsr' - attack-decay-sustain-release type of evolution
%     'curvetype' ('linear')      : type of curve used for non-constant amplitude evolution 
%                                   other options: 'exponential', 'hann-window'
%
%    'attacktime' (0.2)           : attacktime as a percentaje of the duration
%                                   can't be greater than 0.9 in attack-decay mode, or greater than 0.4 in 'adsr' mode
%
% output:
% 	d   - audio signal, sum of harmonics
%	t   - vector of sample time instants
%	f   - matrix where each row is the frequency evolution of each harmonic
%	a   - matrix where each row is the amplitude evolution of each harmonic
%
%
% examples:
%
% [d, t, f, a] = harmonicsource(8000, 1, 4, [500, 800],'f0type', 'linear', 'amplitudes', [0.5, 0.25, 0.20, 0.05])
% four harmonics, linear evolution from 500 to 800 Hz, amplitudes set by user, constant amplitudes, no messages
%
% [d, t, f, a] = harmonicsource(8000, 1, 6, 500, 'f0type', 'constant', 'inharmonicity', 0.001)
% six harmonics, constant f0 of 500 Hz, inharmonicity factor of 0.001, constant amplitudes, no messages
%
% [d, t, f, a] = harmonicsource(8000, 1, 6, 500, 'inharmonicity', 0.001, 'amptype', 'attack-decay', 'attacktime', 0.1)
% same situation of the last example but amplitudes evolve in an attack-decay fashion with linear curves and 10% attack
%
% [d, t, f, a] = harmonicsource(8000, 1, 8, [500, 200],'f0type', 'quadratic', 'verbose', 1)
% eight harmonics, quadratic evolution from 500 to 100 Hz, constant amplitudes, some messages are displayed (e.g f0type)
%
% [d, t, f, a] = harmonicsource(8000, 1, 8, [500, 200],'f0type', 'quadratic', 'amptype', 'adsr', 'curvetype', 'hann-window')
% same situation of the last example but amplitudes evolve in an adsr fashion using a hann-window type of curve
%
% [d, t, f, a] = harmonicsource(8000, 1, 7, [500, 6],'f0type', 'sinusoidal', 'verbose', 2);
% seven harmonics, sinusoidal evolution of 6 Hz centered at 500 Hz, spectrogram and waveform are plotted
%
% Notice that if aliasing occurs a warning is displayed (if verbose mode >= 1)

if nargin < 1;    fs = 8000; end
if nargin < 2;   dur = 1;    end
if nargin < 3; nharm = 4;    end
if nargin < 4;    f0 = 100;  end

% parse out the optional arguments
[verbose, f0type, amplitudes, inharmonicity, amptype, curvetype, attacktime] = ... 
process_options(varargin, 'verbose', 0, 'f0type', 'constant', 'amplitudes', ones(nharm,1), ... 
'inharmonicity', 0, 'amptype', 'constant', 'curvetype', 'linear', 'attacktime', 0.2);

% time
t = 0:1/fs:dur;%t = t(1:end-1);
% fundamental frequency evolution
[f0t, phi] = f0vector(f0, f0type, t, verbose);
% inharmonic factors for each harmonic
inharmonFactors = sqrt(1 + ((1:1:nharm).^2 - 1) * inharmonicity);
% instantaneous frequency evolution of each harmonic
f = ((1:1:nharm) .* inharmonFactors)' * f0t;
% phase evolution of each harmonic
phis = ((1:1:nharm) .* inharmonFactors)' * phi;
% amplitude evolution function
ampt = ampvector(amptype, t, curvetype, attacktime, verbose);
% amplitudes of each harmonic
if length(amplitudes) < nharm
	amplitudes = ones(nharm,1);
	if verbose >= 1; disp('WARNING: not enough amplitude values, forced them all to unity'); end;
end
amps(:,1) = amplitudes / sum(amplitudes);
% amplitude evolution of each harmonic
a = amps * ampt;

% waveform obtained as the sum of harmonics
if (nharm > 1)
	d = sum(a .* sin( 2 * pi * phis));
else
	d = a(1) * sin( 2 * pi * phis);
end

% aliasing warning
if verbose >= 1; if max(max(nharm*f0t)) > fs/2; disp('WARNING: aliasing!'); end; end;

% plots
if (verbose >= 2)
	figure; 
	subplot(3,1,1:2)
	t_window = 0.025;
	window = floor(t_window * fs);
	noverlap = floor(window / 2);
	nfft = window * 2;
	[S,F,T] = spectrogram(d,window,noverlap,nfft,fs);
	imagesc(T,F,10*log10(abs(S)+eps))
	colormap('gray');set(gca,'YDir','normal')
	title('Spectrogram');xlabel('Time (s)');ylabel('Frequency (Hz)');
	subplot(3,1,3)
	plot(t,d,'k');
	title('Waveform');xlabel('Time (s)');ylabel('Amplitude')

	figure; 
	subplot(3,1,1:2)
	plot(t,f','k')
	grid
	ylim([0, fs/2])
	title('Frequency evolution of each harmonic');xlabel('Time (s)');ylabel('Frequency (Hz)');
	subplot(3,1,3)
	plot(t,d,'k');
	title('Waveform');xlabel('Time (s)');ylabel('Amplitude')

	figure; 
	subplot(3,1,1:2)
	plot3(t,f,a,'k')
	grid; ylim([0, fs/2]); view([60 20])
	title('Frequency and amplitude evolution of each harmonic');xlabel('Time (s)');ylabel('Frequency (Hz)');zlabel('Amplitude');
	subplot(3,1,3)
	plot(t,d,'k');
	title('Waveform');xlabel('Time (s)');ylabel('Amplitude')
end

end % end function


function [f0t, phi] = f0vector(f0, f0type, t, verbose)
% instantaneous frequency evolution according to a certain function f0type
	switch f0type
	   case {'constant'}
	      if verbose >= 1; disp('f0 function is constant'); end;
	      phi = f0(1) * t;
	      f0t = f0(1) * ones(1,length(t));
	      
	   case {'linear'}
	      % check if f0 has two items
	      if (length(f0) >= 2)
	      	if verbose >= 1; disp('f0 function is linear'); end;
		phi = (f0(2)-f0(1))/(t(end)-t(1))*((t.^2-t(1)^2)/2+t(1)*(t-t(1)))+f0(1)*(t-t(1));
	      	f0t = (f0(2)-f0(1))/(t(end)-t(1))*(t-t(1))+f0(1);
              else
	      	disp('WARNING: missing f0 value, forced to constant')
	      	[f0t, phi] = f0vector(f0, 'constant', t, verbose);
	      end

	   case 'quadratic'
   	      % check if f0 has two items
        if (length(f0) >= 2)
	      	if verbose >= 1; disp('f0 function is quadratic'); end;
                f0t = (f0(2)-f0(1))/(t(end)-t(1))^2*((t-t(1)).^2)+f0(1);
                phi = (f0(2)-f0(1))/(t(end)-t(1))^2*((t.^3-t(1)^3)/3-t(1)*(t.^2-t(1))+t(1)^2*(t-t(1)))+f0(1)*(t-t(1));
        else
                disp('WARNING: missing f0 value, forced to constant')
                [f0t, phi] = f0vector(f0, 'constant', t, verbose);
        end

	   case 'sinusoidal'
      	      % check if f0 has two items
	      if (length(f0) >= 2)
	      	if verbose >= 1; disp('f0 function is sinusoidal'); end;
                phi = f0(1)*(1-2^(1/(12*2)))/(2*pi*f0(2))*(-cos(2*pi*f0(2)*t)+cos(2*pi*f0(2)*t(1)))+f0(1)*(t-t(1));
                f0t = f0(1)*(1-2^(1/(12*4)))*sin(2*pi*f0(2)*t)+f0(1);
          else
                disp('WARNING: missing f0 value, forced to constant')
                [f0t, phi] = f0vector(f0, 'constant', t, verbose);
	      end
	      
	   otherwise
	      disp('WARNING: unknown f0 function, forced to constant')
	      f0t = f0vector(f0, 'constant', t, verbose)
	end

end


function ampt = ampvector(amptype, t, curvetype, attacktime, verbose)
% instantaneous frequency evolution according to a certain function f0type
maxamp = 1;
% check attacktime
if (attacktime > 0.9); attacktime=0.2; disp('WARNING: attack time larger than 90%, forced to 20\%'); end;
% attack point
k = find(t > attacktime * t(end), 1, 'first');

	switch amptype
	   case {'constant'}
	      if verbose >= 1; disp('amplitude function is constant'); end;
	      ampt = ones(1,length(t));
     
	   case {'attack-only'}
		if verbose >= 1; disp(['amplitude function is ' amptype ' ' curvetype]); end;
		switch curvetype
			case {'linear'}
				ampt = maxamp * 1/(t(end)-t(1)).*(t-t(1));
			case {'exponential'}
				ampt = maxamp * 1/(t(end)-t(1))^2.*(t-t(1)).^2;
			case {'hann-window'}
			   	ampt = maxamp * hann(2*length(t))';
				ampt = ampt(1:end/2);
			otherwise
				disp('WARNING: unknown curve type, forced to linear')
				ampt = ampvector(amptype, t, 'linear', attacktime, verbose)
		end
		
	   case {'decay-only'}
		if verbose >= 1; disp(['amplitude function is ' amptype ' ' curvetype]); end;
		switch curvetype
			case {'linear'}
				ampt = maxamp * -1/(t(end)-t(1)).*(t-t(1)) + maxamp;
			case {'exponential'}
				ampt = maxamp * -1/(t(end)-t(1))^2.*(t-t(1)).^2 + maxamp;
			case {'hann-window'}
			   	ampt = maxamp * hann(2*length(t))';
				ampt = ampt(end/2+1:end);
			otherwise
				disp('WARNING: unknown curve type, forced to linear')
				ampt = ampvector(amptype, t, 'linear', attacktime, verbose)
		end
				
	   case {'attack-decay'}
		if verbose >= 1; disp(['amplitude function is ' amptype ' ' curvetype]); end;
		switch curvetype
			case {'linear'}
				ampt = maxamp * [1/(t(k)-t(1)).*(t(1:k)-t(1)) -1/(t(end)-t(k+1)).*(t(k+1:end)-t(k+1)) + maxamp];
			case {'exponential'}
				ampt = maxamp * [1/(t(k)-t(1))^2.*(t(1:k)-t(1)).^2 -1/(t(end)-t(k+1))^2*(t(k+1:end)-t(k+1)).^2+maxamp];
			case {'hann-window'}
			   	ampt = maxamp * hann(length(t))';
			otherwise
				disp('WARNING: unknown curve type, forced to linear')
				ampt = ampvector(amptype, t, 'linear', attacktime, verbose)
		end

	   case {'adsr'}
		if verbose >= 1; disp(['amplitude function is ' amptype ' ' curvetype]); end;
		if (attacktime >= 0.4)
			attacktime = 0.2;
			k = find(t > attacktime * t(end), 1, 'first');
			disp('WARNING: attack time too large, forced to 20%'); 
		end;
		d = floor(3*k/2);
		s = length(t)-k;
		switch curvetype
			case {'linear'};
				ampt = [maxamp/(t(k)-t(1)).*(t(1:k)-t(1)) (2*maxamp/3-maxamp)/(t(d)-t(k+1)).*(t(k+1:d)-t(k+1))+maxamp 2*maxamp/3*ones(1,length(d+1:s)) -2*maxamp/3/(t(end)-t(s+1)).*(t(s+1:end)-t(s+1))+2*maxamp/3];
			case {'exponential'}
				ampt = [maxamp/(t(k)-t(1))^2.*(t(1:k)-t(1)).^2 (2*maxamp/3-maxamp)/(t(d)-t(k+1))^2.*(t(k+1:d)-t(k+1)).^2+maxamp 2*maxamp/3*ones(1,length(d+1:s)) -2*maxamp/3/(t(end)-t(s+1))^2.*(t(s+1:end)-t(s+1)).^2+2*maxamp/3];
			case {'hann-window'}
				ampt = hann(2*length(t(1:k)))';
				atmp = maxamp/3*hann(2*length(t(k+1:d)))';
				ampt = [maxamp*ampt(1:k) atmp(length(t(k+1:d))+1:end)+2*maxamp/3 2*maxamp/3*ones(1,length(t(d+1:s))) 2*maxamp/3*ampt(k+1:end)];
			otherwise
				disp('WARNING: unknown curve type, forced to linear')
				ampt = ampvector(amptype, t, 'linear', attacktime, verbose)
		end
   
	   otherwise
	      disp('WARNING: unknown amptype, forced to constant')
	      ampt = ampvector('constant', t, curvetype, attacktime, verbose)
	end

end

% PROCESS_OPTIONS - Processes options passed to a Matlab function.
%                   This function provides a simple means of
%                   parsing attribute-value options.  Each option is
%                   named by a unique string and is given a default
%                   value.
%
% Usage:  [var1, var2, ..., varn[, unused]] = ...
%           process_options(args, ...
%                           str1, def1, str2, def2, ..., strn, defn)
%
% Arguments:   
%            args            - a cell array of input arguments, such
%                              as that provided by VARARGIN.  Its contents
%                              should alternate between strings and
%                              values.
%            str1, ..., strn - Strings that are associated with a 
%                              particular variable
%            def1, ..., defn - Default values returned if no option
%                              is supplied
%
% Returns:
%            var1, ..., varn - values to be assigned to variables
%            unused          - an optional cell array of those 
%                              string-value pairs that were unused;
%                              if this is not supplied, then a
%                              warning will be issued for each
%                              option in args that lacked a match.
%
% Examples:
%
% Suppose we wish to define a Matlab function 'func' that has
% required parameters x and y, and optional arguments 'u' and 'v'.
% With the definition
%
%   function y = func(x, y, varargin)
%
%     [u, v] = process_options(varargin, 'u', 0, 'v', 1);
%
% calling func(0, 1, 'v', 2) will assign 0 to x, 1 to y, 0 to u, and 2
% to v.  The parameter names are insensitive to case; calling 
% func(0, 1, 'V', 2) has the same effect.  The function call
% 
%   func(0, 1, 'u', 5, 'z', 2);
%
% will result in u having the value 5 and v having value 1, but
% will issue a warning that the 'z' option has not been used.  On
% the other hand, if func is defined as
%
%   function y = func(x, y, varargin)
%
%     [u, v, unused_args] = process_options(varargin, 'u', 0, 'v', 1);
%
% then the call func(0, 1, 'u', 5, 'z', 2) will yield no warning,
% and unused_args will have the value {'z', 2}.  This behaviour is
% useful for functions with options that invoke other functions
% with options; all options can be passed to the outer function and
% its unprocessed arguments can be passed to the inner function.

% Copyright (C) 2002 Mark A. Paskin
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
% General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
% USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = process_options(args, varargin)

% Check the number of input arguments
n = length(varargin);
if (mod(n, 2))
  error('Each option must be a string/value pair.');
end

% Check the number of supplied output arguments
if (nargout < (n / 2))
  error('Insufficient number of output arguments given');
elseif (nargout == (n / 2))
  warn = 1;
  nout = n / 2;
else
  warn = 0;
  nout = n / 2 + 1;
end

% Set outputs to be defaults
varargout = cell(1, nout);
for i=2:2:n
  varargout{i/2} = varargin{i};
end

% Now process all arguments
nunused = 0;
for i=1:2:length(args)
  found = 0;
  for j=1:2:n
    if strcmpi(args{i}, varargin{j})
      varargout{(j + 1)/2} = args{i + 1};
      found = 1;
      break;
    end
  end
  if (~found)
    if (warn)
      warning(sprintf('Option ''%s'' not used.', args{i}));
      args{i}
    else
      nunused = nunused + 1;
      unused{2 * nunused - 1} = args{i};
      unused{2 * nunused} = args{i + 1};
    end
  end
end

% Assign the unused arguments
if (~warn)
  if (nunused)
    varargout{nout} = unused;
  else
    varargout{nout} = cell(0);
  end
end
 
end