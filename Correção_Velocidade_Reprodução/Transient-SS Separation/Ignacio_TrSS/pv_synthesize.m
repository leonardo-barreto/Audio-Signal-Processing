function X = pv_synthesize(M, P, win, synt_hop, an_hop)


[ num_bins, num_frames ] = size(P);

delta_phi=zeros(num_bins, num_frames-1);

PF=zeros(num_bins, num_frames);

%window=hanningz(win); % tapering window
window=hanning(win); % tapering window
%

% phase unwrapping

%
two_pi=2*pi;
omega = two_pi*an_hop*[0:num_bins-1]'/num_bins;
for idx=2 : num_frames
    
    ddx = idx-1;
    
    delta_phi(:,ddx) =princarg(P(:,idx)-P(:,ddx)-omega);
    
    phase_inc(:,ddx)=(omega+delta_phi(:,ddx))/an_hop;
end
%
% now prepare a matrix of complex numbers which
% recombine modulo and phases to be able to feed
% the ifft algorithms
%
% the idea here is to use the values of the previous
% phases, and calculate the current phases computing
% the current phase difference multiplied by the
% current hop size
%


PF(:,1)=P(:,1); % the initial phase is the same
for idx = 2:num_frames
    
    ddx = idx-1;
    PF(:,idx)=PF(:,ddx)+synt_hop*phase_inc(:,ddx);
end;
Z=M.*exp(i*PF);
%
% perform inverse windowing and overlap-adding
% of the resulting ifft frames
%
textprogressbar('Sinthetizing: ')
X = zeros((num_frames*synt_hop)+win, 1);
curstart = 1;
for idx = 1:num_frames
    curend = curstart + win - 1;
    RIfft = fftshift(real(ifft(Z(:,idx))));
    X([curstart:curend])=X([curstart:curend])+RIfft.*window;
    curstart = curstart + synt_hop;
%fprintf('\b'); % delete previous counter display
textprogressbar(idx/num_frames*100);

end

k=sum(hanning(win) .* window)/synt_hop;
X=X/k;
textprogressbar('Done.')

end
