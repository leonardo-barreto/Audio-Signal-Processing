function [C, v1, v2] = structure_tensor_v2(varargin)
%STRUCTURE_TENSOR_V2 computes the structure tensor of the input matrix
%using the X > 0 criterion to the anisotropy computation.
%   Detailed explanation goes here

X_structure_tensor = varargin{1};
Ngauss = varargin{2};
sigma_m = varargin{3};
sigma_k = varargin{4};

% Computing the partial derivatives
[X_m, X_k] = imgradientxy(X_structure_tensor);

% Computing T11, T21, T12 and 22
G = assym_gauss2D(Ngauss, [sigma_m sigma_k], 1); % Using symmetric

X_mm = conv2(X_m.*X_m, G, 'same');
X_mk = conv2(X_k.*X_m, G, 'same'); X_km = X_mk;
X_kk = conv2(X_k.*X_k, G, 'same');

% Computing angles and anisotropy measure
C = zeros(size(X_structure_tensor));
v1 = zeros(size(X_structure_tensor));
v2 = zeros(size(X_structure_tensor));

for m = 1:size(X_structure_tensor,2) % time frames
    for k = 1:size(X_structure_tensor ,1) % frequency bins
        % Anisotropy measure
        if X_structure_tensor(k,m) > 0
            T = [X_mm(k,m) X_km(k,m); X_mk(k,m) X_kk(k,m)];

            [eig_vec, eig_val] = eig(T);

            % Angle of orientation
            v1(k,m) = eig_vec(1,1);
            v2(k,m) = eig_vec(2,1);
        
            C(k,m) = ((eig_val(4) - eig_val(1)) / (eig_val(1) + eig_val(4)))^2;
        end
    end
end

% Smoothing C transitions
G_C = assym_gauss2D(5, [1 1], 1); % Using symmetric
C = conv2(C, G_C, 'same');

end

function G = assym_gauss2D(N, sigma, assym_factor)
% N should be odd
N = ceil(N);

if ~mod(N,2)
    N = N - 1;
end

x  = -(N-1)/2:(N-1)/2;

gm = 1/sqrt(2*pi*sigma(1)^2)*exp(-x.^2/(2*sigma(1)^2)); gm = gm(:);
gm_2 = 1/sqrt(2*pi*sigma(1)^2)*exp(-x.^2/(2*(sigma(1)*assym_factor)^2)); gm_2 = gm_2(:);
gm = [gm_2(1:(end-1)/2); gm((end-1)/2+1:end)];

gk = 1/sqrt(2*pi*sigma(2)^2)*exp(-x.^2/(2*sigma(2)^2)); gk = gk(:);

G = gk*(gm');

end
