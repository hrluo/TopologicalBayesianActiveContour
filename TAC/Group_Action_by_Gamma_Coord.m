function fn = Group_Action_by_Gamma_Coord(f,gamma)
% Inputs:
% f = curve stored as n*T matrix, NOT SRVF
% gamma = re-parameterization function over [0,1]

% Output:
% fn = re-parameterization of f by gamma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(f);

% Re-parameterize f
for j=1:n
    fn(j,:) = interp1(linspace(0,1,T),f(j,:),gamma);
end

return;