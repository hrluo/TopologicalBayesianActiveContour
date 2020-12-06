function gamI = invertGamma(gam)
% Input:
% gam = gamma values at certain points

% Output:
% gamI = inverse of gam (necessary since output of DynamicProgrammingQ.c
% gives the inverse re-parameterization function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = length(gam);

% Compute inverse of gam assumed to be uniformly spaced on [0,1] with N
% steps
x = [1:N]/N;
gamI = interp1(gam,x,x,'linear');
