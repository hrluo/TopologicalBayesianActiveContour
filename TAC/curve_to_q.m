function q = curve_to_q(p)
% Input:
% n x N matrix with columns containing coordinates of the curve p

% Output:
% n x N matrix representing the unit-norm SRVF of curve p
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,N] = size(p);

% Compute time-derivative of p
for i = 1:n
    v(i,:) = gradient(p(i,:),1/N);
end

% Compute SRVF
for i = 1:N
    L(i) = sqrt(norm(v(:,i),'fro'));
    if L(i) > 0.00001
        q(:,i) = v(:,i)/L(i);
    else
        q(:,i) = 0*ones(n,1);
    end
end

% Re-scale SRVF to be unit norm
q = q/sqrt(InnerProd_Q(q,q));