function delG = Form_Basis_Normal_A(q)
% Input:
% q = SRVF of curve for which we want an orthonormal basis of vectors for
% its normal space, with respect to inner product

% Output:
% delG = basis vectors (as a cell array)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(q);
e = eye(n);
for i = 1:n
    Ev(:,:,i) = repmat(e(:,i),1,T);
end

for t = 1:T
    qnorm(t) = norm(q(:,t));
end

% Compute basis vectors
for i = 1:n
    tmp1 = repmat(q(i,:)./qnorm,n,1);
    tmp2 = repmat(qnorm,n,1);
    delG{i} = tmp1.*q + tmp2.*Ev(:,:,i);    
end

return;