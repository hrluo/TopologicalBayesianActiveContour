function K = FindElasticCovariance(q_bar,q,reparamFlag)
% Inputs:
% q_bar = 2 x N matrix of estimated mean shape SRVF
% q = 2 x N x M matrix of SRVF coordinates, where N = number of
% discretization points and M = number of SRVFs in sample
% reparamFlag = 1 if want to find optimal re-parameterization

% Output:
% K = covariance matrix on tangent space at q_bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b,n] = size(q);

% Default to searching for optimal re-parameterization if unspecified
if ~exist('reparamFlag','var') || isempty(reparamFlag), reparamFlag = 1; end

for i=1:n
    % Project SRVFs into tangent space at q_bar
    tmp = ElasticShootingVector(q_bar,q(:,:,i),reparamFlag);
    
    % Stack 2 x N tangent vector as 2N column vector
    VV(i,1:b) = tmp(1,:);
    VV(i,b+1:2*b) = tmp(2,:);
end

% Compute covariance matrix of vectorized tangent vectors on tangent space
% at q_bar
K = cov(VV);