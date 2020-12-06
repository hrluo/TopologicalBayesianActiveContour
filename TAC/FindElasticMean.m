function [beta_bar,q_bar] = FindElasticMean(Data,figs)
% Inputs:
% Data = 2 x N x M vector-valued curve coordinates (NOT SRVF)
% figs = 1 if want to display the mean shape and energy plot (0 otherwise)

% Outputs:
% beta_bar = estimated mean shape
% q_bar = estimated mean shape SRVF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Niter = 20;                 % number of iterations to run gradient-descent
[d,N,n] = size(Data);

% Compute SRVFs of all curves
for i=1:n
    X = ReSampleCurve(Data(:,:,i),N);
    q(:,:,i) = curve_to_q(X);
end

% Initialize gradient descent algorithm at pointwise mean
del = 0.5;                  % step size for tangent vector update
q_bar = sum(q,3)/n;
q_bar = q_bar/sqrt(InnerProd_Q(q_bar,q_bar));
q_bar = ProjectC(q_bar);

% Iterate to find Karcher mean
for iter=1:Niter
    vm = 0;    
    for i=1:n
        [iter i]
        [v,~,~] = ElasticShootingVector(q_bar,q(:,:,i),1);
        vm = vm+v;
    end
       vm = vm/n;
       
       E(iter) = norm(vm,'fro');
       q_bar = ElasticShooting(q_bar,del*vm);
       
       iter = iter+1;
       E
    
end

% Convert mean SRVF to curve representation
beta_bar = q_to_curve(q_bar);

% Display energy and mean shape plots, if desired
if figs==1
    figure(101);
    plot(E);
    
    figure(21); clf;
    plot(beta_bar(1,:),beta_bar(2,:),'LineWidth',3);
    axis equal; axis off;
end