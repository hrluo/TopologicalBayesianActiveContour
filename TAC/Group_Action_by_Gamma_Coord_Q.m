function qn = Group_Action_by_Gamma_Coord_Q(q,gamma)
% Inputs:
% q = SRVF stored as n*T matrix
% gamma = reparameterization function over [0,1]

% Output:
% qn = re-parameterized SRVF of q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(q);

% Compute time-derivative of gamma
gamdot = gradient(gamma,1/T);

% Re-parameterize q by gamma
for j=1:n
    qn(j,:) = interp1(linspace(0,1,T),q(j,:),gamma);
    qn(j,:)=qn(j,:).*sqrt(gamdot);
end

return;