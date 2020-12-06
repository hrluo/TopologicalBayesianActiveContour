function val = InnerProd_Q(q1,q2)
% Inputs:
% q1, q2 = n x T matrix of SRVFs for each of two curves

% Output:
% val = L2 inner product of two curves q1, q2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(q1);

% Compute L2 inner product
val = trapz(linspace(0,1,T),sum(q1.*q2));