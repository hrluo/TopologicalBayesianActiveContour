function p = q_to_curve(q)
% Input:
% n x N matrix with columns containing coordinates of the curve p

% Output:
% n x N matrix representing the unit-norm SRVF of curve p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,N] = size(q);

% Compute Euclidean norm for each coordinate
for i=1:N
    qnorm(i) = norm(q(:,i),'fro');
end

% Map SRVF q to curve p
for i=1:n
    p(i,:) = [cumtrapz(q(i,:).*qnorm )/N] ;
end

return;