function [delP,beta_shift] = PriorUpdate(beta,q_bar,A,figs)
% Inputs:
% beta = current curve contour
% q_bar = SRVF of average of training curves
% A = matrix update computed based on svd of covariance matrix
% figs = 1 if want to see prior update at each step

% Output:
% delP = prior update
% beta_shift = optimal shift of "starting point" of beta to best match
% shape given by q_bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,N] = size(beta);

% Calculate and save length and centroid of current curve
% Length
for i=1:d
    gradb(i,:) = gradient(beta(i,:),1/(N-1));
end
L = trapz(linspace(0,1,N),vecnorm(gradb));

% Centroid
cen = mean(beta,2);

% Convert curve to SRVF
q = curve_to_q(beta);

% Obtain shooting vector via inverse exponential map
[wv,~,qa,gamI,O,shiftbest] = ElasticShootingVector(q_bar,q,1);
w = [wv(1,:),wv(2,:)]';

% Manually shift "starting point" of beta and apply optimal gamma
beta_shift = ShiftF(beta,shiftbest);
beta_shift_rep = Group_Action_by_Gamma_Coord(beta_shift,gamI);
q_shift_rep = curve_to_q(beta_shift_rep);

% Compute gradient
gradvv = A*w;
gradv(1,:) = gradvv(1:N);
gradv(2,:) = gradvv((N+1):(2*N));

% Parallel transport gradient vector to tangent space at qa
v = Parallel_Transport_C(gradv,q_bar,qa);

% Shoot current contour in direction of negative gradient vector
eps = 0.3;
q_new = ElasticShooting(q_shift_rep,-eps*v);

% Apply inverse of gamma to update to be matched with arc-length
% parameterization of current contour
gam_inv = invertGamma(gamI);
gam_inv = (gam_inv-gam_inv(1))/(gam_inv(end)-gam_inv(1));
beta_new = q_to_curve(q_new);
beta_new = Group_Action_by_Gamma_Coord(beta_new,gam_inv);

% Re-scale and match centroid for updated beta
beta_new = L*beta_new;
beta_new = beta_new-mean(beta_new,2)+cen;

% Update term
delP = (beta_shift-beta_new)/eps;

% Test out
if figs==1
    figure(4)
    clf
    axis equal
    hold on
    plot(beta_shift(1,:),beta_shift(2,:),'b')
    plot(beta_new(1,:),beta_new(2,:),'r')
    for i=1:N
        j = i;
        plot([beta_shift(1,j) beta_new(1,j)],[beta_shift(2,j) beta_new(2,j)],'k')
    end
    plot(beta_shift(1,1),beta_shift(2,1),'b*')
    plot(beta_new(1,1),beta_new(2,1),'r*')
end