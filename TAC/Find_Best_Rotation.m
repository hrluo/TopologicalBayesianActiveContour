function [q2new,R] = Find_Best_Rotation(q1,q2)
% Inputs:
% q1, q2 = SRVFs of two curves

% Outputs:
% q2new = optimally-rotated version of q2 to q1
% R = optimal rotation matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assumes the starting points are fixed

% Basics
[n,~] = size(q1);

% Find optimal rotation matrix
A = q1*q2';
[U,~,V] = svd(A);
if det(A) > 0
    S = eye(n);
else
    S = eye(n);
    S(:,end) = -S(:,end);   % change sign of last column (see Dryden and Mardia (2016))
end
R = U*S*V';

% Apply optimal rotation matrix to q2
q2new = R*q2;