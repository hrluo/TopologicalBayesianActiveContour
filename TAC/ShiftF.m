function pn = ShiftF(p,tau)
% Inputs:
% p = curve matrix stored in 2*T (only use for CLOSED curves)
% tau = the unit of space we want to shift

% Output:
% pn = the curve matrix of p shifted by tau units
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(p);
if(tau == 0)
    pn = p;
    return;
end

% Permute order that coordinates are presented in p according to tau
if tau > 0
    pn(:,1:T-tau) = p(:,tau+1:T);
    pn(:,T-tau+1:T) = p(:,1:tau);
    pn(:,T-tau) = [];       % remove repeated point
    pn(:,end+1) = pn(:,1);
    return;
else
    t = abs(tau)+1;
    pn(:,1:T-t+1) = p(:,t:T);
    pn(:,T-t+2:T) = p(:,1:t-1);
    pn(:,T-t+1) = [];       % remove repeated point
    pn(:,end+1) = pn(:,1);
    return;
end