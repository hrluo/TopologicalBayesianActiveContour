function q2n = ElasticShooting(q1,v)
% Inputs:
% q1 = SRVF of a curve (pole of tangent space)
% v = tangent vector on tangent space at q_1

% Output:
% q2n = SRVF of resulting curve after shooting v on tangent space at q1
% (i.e., exponential map q2n = exp_{q1}(v))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute norm of shooting vector
d = sqrt(InnerProd_Q(v,v));
eps = 0.1;          % step size for computing exponential map

% Find resulting SRVF after moving in the direction of tangent vector v on
% tangent space at q1
if d < 0.00001      % tangent vector is zero vector
    q2n = q1;
else
    % Compute exponential map in a step-wise manner
    q2{1} = cos(eps*d)*q1 + (sin(eps*d)/d)*v;
    q2{1} = ProjectC(q2{1});
    v = Parallel_Transport_C(v,q1,q2{1});
    d = sqrt(InnerProd_Q(v,v));
    for j=2:10
        q2{j} = cos(eps*d)*q2{j-1} + (sin(eps*d)/d)*v;
        q2{j} = ProjectC(q2{j});
        v = Parallel_Transport_C(v,q2{j-1},q2{j});
        d = sqrt(InnerProd_Q(v,v));
    end
    q2n = q2{10};   % final SRVF
end