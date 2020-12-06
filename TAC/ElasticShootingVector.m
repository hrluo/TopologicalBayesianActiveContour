function [v,d,q2n,gamI,R,shiftbest] = ElasticShootingVector(q1,q2,reparamFlag)
% Inputs:
% q1, q2 = SRVFs of two curves
% reparamFlag = 1 if want to find optimal re-parameterization

% Outputs:
% v = tangent vector representing q2n on tangent space at q1
% (i.e., inverse-exponential map v = exp^{-1}_{q1}(q2n))
% d = geodesic distance between q1 and q2n
% q2n = aligned SRVF of second curve to first curve's SRVF
% gamI = optimal re-parameterization function
% R = optimal rotation matrix
% shiftbest = the permutation-based shift of coordinates of q2 which match
% q1 best (i.e., dictates "best" starting point of q2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find optimal rotation and re-parameterization (if reparamFlag = 1)
[q2n,gamI,R,shiftbest] = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag);

q2n = q2n/sqrt(InnerProd_Q(q2n,q2n));

% Compute geodesic distance
d = acos(InnerProd_Q(q1,q2n));

% Compute tangent vector on tangent space at q_1
if d < 0.0001
    v = zeros(size(q1));
else
    v = (d/sin(d))*(q2n - cos(d)*q1);
end