function [E,d] = PriorEnergy(beta,q_bar,Um,Sm_inv,delta)
% Inputs:
% beta = current curve contour
% q_bar = SRVF of average of training curves
% Um = matrix with columns which are eigenvectors obtained from SVD
% of covariance matrix found in FindElasticCovariance.m
% Sm = diagonal matrix of eigenvalues obtained from SVD of covariance
% matrix found in FindElasticCovariance.m
% delta = parameter (typically set to be smaller than the smallest
% eigenvalue found in Sm)

% Outputs:
% E = prior energy term
% d = elastic shape distance between beta and curve with SRVF q_bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert curve to SRVF (which lies on unit Hilbert sphere)
q = curve_to_q(beta);

% Obtain shooting vector via inverse-exponential map
[v,d,~] = ElasticShootingVector(q_bar,q,1);

% Vectorize the shooting vector to match dimension of Um, Sm
w = [v(1,:),v(2,:)]';

% Compute prior energy term
E = 0.5*w'*(Um*Sm_inv*Um')*w + (1/(2*delta^2))*norm(w-Um*Um'*w)^2;