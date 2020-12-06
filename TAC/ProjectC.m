function qnew = ProjectC(q)
% Input:
% q = SRVF stored as n*T matrix

% Output:
% qnew = projection of SRVF q into space of closed curves (even if q is
% open to begin with)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(q);

% Set spacing
if(n == 2)
    dt = 0.35;
end
if(n == 3)
    dt = 0.2;
end
epsilon = 1e-6;

e = eye(n);
iter = 1;
res = ones(1,n);
J = eye(n);

s = linspace(0,1,T);
qnew = q;

C = [];
while (norm(res,2) > epsilon)
    if(iter > 300)
        iter
        break;
    end
    
    % Compute Jacobian
    for i=1:n
        for j=1:n
            J(i,j) = 3 * trapz(s,qnew(i,:) .*qnew(j,:) );
        end
    end
    

    for i=1:T
        qnorm(i) = norm(qnew(:,i),2);
    end

    % Compute the residue
    for i = 1:n
        G(i) = trapz(s,qnew(i,:).*qnorm);
    end
    res = -G;

    if(norm(res) < epsilon)
        break;
    end

    x = inv(J)*res';
    C(iter) = norm(res,2);

    delG = Form_Basis_Normal_A(qnew);
    temp = 0;
    for i = 1:n
        temp = temp + x(i)*delG{i}*dt;
    end
    qnew = qnew + temp;

    iter = iter + 1;
end