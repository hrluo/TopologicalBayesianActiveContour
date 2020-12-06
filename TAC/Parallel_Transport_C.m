function w_new = Parallel_Transport_C(w,q1,q2)
% Inputs:
% w = tangent vector on tangent space at q1
% q1 = SRVF which is the pole of the tangent space w lies in
% q2 = SRVF which is the pole of the tangent space we are parallel
% transporting w to

% Output:
% w_new = tangent vector on tangent space at q2 obtained by parallel
% transporting w
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lw = sqrt(InnerProd_Q(w,w));        % norm of tangent vector w

if(lw < 0.0001)
    w_new = w;                      % no transportation needed
else
    w_new = Project_Tangent(w,q2);
    w_new = w_new*lw/sqrt(InnerProd_Q(w_new,w_new));
end

return;
