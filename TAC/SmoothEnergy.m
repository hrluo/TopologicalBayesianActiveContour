function E = SmoothEnergy(beta)
% Input:
% beta = curve

% Output:
% E = smoothing energy term
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting
[d,N] = size(beta);

% Compute smoothing energy functional
for i=1:d
    v(i,:) = gradient(beta(i,:),1/(N-1));
end
E = trapz(linspace(0,1,N),vecnorm(v,2,1));