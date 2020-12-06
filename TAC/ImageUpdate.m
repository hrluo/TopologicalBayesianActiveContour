function sc_v = ImageUpdate(beta,llikI)
% Inputs: 
% beta = current curve contour
% llikI = log-likelihood image

% Output:
% sc_v = univariate function of negative log-likelihoods along current
% contour (will be multiplied by outward unit normal vector for update)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
N = size(beta,2);

% Update
for j=1:N
    sc_v(j) = -llikI(ceil(beta(1,j)),ceil(beta(2,j)));
end