function nrm_v = OutwardUnitNormal(beta)
% Input:
% beta = 2 x N matrix of coordinates for curve

% Output:
% nrm_v = outward normal vector to beta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(beta,2);   % number of sampling points along curve

% Compute gradients in each direction
t(1,:) = gradient(beta(1,1:(N-1)),1/(N-2));
t(2,:) = gradient(beta(2,1:(N-1)),1/(N-2));
t(:,N) = t(:,1);

% Outward normal vector field
nrm_v = [-t(2,:);t(1,:)]./(vecnorm(t));