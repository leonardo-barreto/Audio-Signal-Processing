function [vout] = TimeVaryingResample (vin, facvec, offvec)
% tic
sincTail = 100;

facveclen = length(facvec);
vinlen = length(vin);
vout = [];
if(facveclen ~= length(offvec)),
  error('ERROR: resampling factors vector length and offset vector length must be the same.')
elseif (offvec(end) > vinlen),
  error('Erro! Last offset greater than input length.')
end
counter = 1;    
dev = 0;    % deviation caused by number of samples approximation 

vout(1:(offvec(1)- 1)) = vin(1:(offvec(1)- 1))';
vout = vout';
vinTemp = [vin;zeros(sincTail,1)];

while counter < facveclen
    ta = 1/facvec(counter);         % initial resampling factor
    tb = 1/facvec(counter + 1);     % final resampling factor
    offset = offvec(counter);           % time corresponding to initial resampling factor
    nextOffset = offvec(counter + 1);   % time corresponding to final resampling factor
    samplesIn = nextOffset - offset;    % number of samples to resample
    
    a = ta+tb; b = 3*ta-tb+2-2*samplesIn; c = 2-2*samplesIn;   % 2nd order equation parameters
    samplesOut = floor(max(roots([a b c])))+1;   % number of samples after resampling
    k = (tb-ta)/(samplesOut);                 % time-varying resampling step
        
    step = cumsum([0 (0:k:k*(samplesOut-2))]) + (0:ta:(samplesOut-1)*ta);
    positions = offset + step - dev;
    centers = round(positions)';                                 % rounded values for each sinc center
    e = positions' - centers;                  % rounding error

    mat1 = repmat((-sincTail:sincTail) , samplesOut,1);     % nyquist filter matrix
    mat2 = repmat(e,1,2*sincTail+1);                        % matrix for error compensation
    matSinc = sinc(pi*(mat1-mat2));                         % sinc matrix with error compensation
    
    mat3 = repmat(centers,1,2*sincTail+1);  % shifting to each output sample
    
    block = dot(matSinc,vinTemp(mat1 + mat3),2);        % vout is just a simple dot product row wise
    
%     procvec = vinTemp(offset + (-sincTail:(samplesIn + sincTail -1)));
%     a = SlopeResample(procvec,ta,k,samplesOut,dev,sincTail);
    vout = [vout; block];
    dev = dev + samplesIn - samplesOut*ta - samplesOut*(samplesOut-1)*k/2;
    counter = counter + 1;
end
offset = offvec(counter);
step = (0:1:(vinlen-offset-1));
positions = offset + step - dev;
centers = round(positions)';                                 % rounded values for each sinc center
e = positions' - centers;                  % rounding error

mat1 = repmat((-sincTail:sincTail) , vinlen-offset,1);     % nyquist filter matrix
mat2 = repmat(e,1,2*sincTail+1);                        % matrix for error compensation
matSinc = sinc(pi*(mat1-mat2));                         % sinc matrix with error compensation

mat3 = repmat(centers,1,2*sincTail+1);  % shifting to each output sample
finalBlock = dot(matSinc,vinTemp(mat1 + mat3),2);        % vout is just a simple dot product row wise
vout = [vout; finalBlock];
% toc
end

% while counter < facveclen
%     ta = 1/facvec(counter);         % initial resampling factor
%     tb = 1/facvec(counter + 1);     % final resampling factor
%     offset = offvec(counter);           % time corresponding to initial resampling factor
%     nextOffset = offvec(counter + 1);   % time corresponding to final resampling factor
%     samplesIn = nextOffset - offset;    % number of samples to resample
%     procvec = vinTemp(offset + (-sincTail:(samplesIn + sincTail -1)));
% 
% %     a = ta+tb; b = 3*ta-tb+2-2*samplesIn; c = 2-2*samplesIn;   % 2nd order equation parameters
% %     samplesOut = floor(max(roots([a b c])))+1;   % number of samples after resampling
% %     k = (tb-ta)/(samplesOut);                 % time-varying resampling step
%     
% %     if counter > 500
% %         b = 3;
% %     end
% 
%     samplesOut = floor((2*samplesIn - 2)/(ta+tb)) + 1;
%     k = (tb-ta)/(samplesOut - 2);
%     
%     a = SlopeResample(procvec,ta,k,samplesOut,dev,sincTail);
%     vout = [vout; a];
% %     dev = dev + samplesIn - samplesOut*ta - samplesOut*(samplesOut-1)*k/2;
%     dev = dev + samplesIn - (samplesOut)*ta - (samplesOut+1)*(samplesOut-2)*k/2;
%     counter = counter + 1;
% end
% procvec = vinTemp((offvec(counter)-sincTail):end);
% vout = [vout; SlopeResample(procvec,1,0,vinlen-offvec(counter),dev,sincTail)];
% % toc
% end