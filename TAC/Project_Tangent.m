function fnew = Project_Tangent(f,q)
% Inputs:
% f = tangent vector to be projected into tangent space at q
% q = SRVF represent pole of tangent space

% Output:
% fnew = projection of f into tangent space at q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(q);

% Project w in T_q({\cal B}), ie the unit sphere
w = f - InnerProd_Q(f,q)*q;
e = eye(n);

% Form the basis for the Normal space of {\cal A}
g = Form_Basis_Normal_A(q);

% Refer to the function Gram_Schmidt for the parameters
Evorth = Gram_Schmidt(g,'InnerProd_Q');

% Unpack Evorth structure
for i=1:n
    Ev(:,:,i) = Evorth{i};
end

for i=1:n
    w = w - InnerProd_Q(w,Ev(:,:,i))*Ev(:,:,i);
end

fnew=w;

return;