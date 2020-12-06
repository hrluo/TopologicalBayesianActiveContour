function [q2best,gamIbest,Rbest,shiftbest] = Find_Rotation_and_Seed_unique(q1,q2,reparamFlag)
% Inputs:
% q1, q2 = SRVFs of two curves
% reparamFlag = 1 if want to find optimal re-parameterization

% Outputs:
% q2best = aligned SRVF of second curve to first curve's SRVF
% gamIbest = optimal re-parameterization function
% Rbest = optimal rotation matrix
% shiftbest = the permutation-based shift of coordinates of q2 which match
% q1 best (i.e., dictates "best" starting point of q2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,T] = size(q1);

% Default to searching for optimal re-parameterization if unspecified
if ~exist('reparamFlag','var') || isempty(reparamFlag), reparamFlag = 1; end

% Other settings
scl = 1;    % how often to check for best "starting point" of q2 (i.e.,
            % setting scl = 2 checks every 2nd point of q2 to determine
            % if it is the best "starting point")
minE = 1000;    % initialize value for minimal energy

% Find optimal rotation, re-parameterization at each "starting point" and
% select one which minimizes the elastic distance
for ctr = 0:floor((T-1)/scl)
    % Permute order of coordinates in q2
    q2n = ShiftF(q2,scl*ctr);
    
    % Find optimal rotation
    [q2n,R] = Find_Best_Rotation(q1,q2n);
    
    % Find optimal re-parameterization
    if(reparamFlag)    
        if norm(q1-q2n,'fro') > 0.0001
            gam = DynamicProgrammingQ(q1,q2n,0,0);
            gamI = invertGamma(gam);
            gamI = (gamI-gamI(1))/(gamI(end)-gamI(1));
            q2new = Group_Action_by_Gamma_Coord_Q(q2n,gamI);
            q2new = ProjectC(q2new);
            q2new = q2new/sqrt(InnerProd_Q(q2new,q2new));
        else
            q2new = q2n;    % q1 is approximately equal to q2n
        end
    else
        q2new  = q2n;
    end

    % Find optimal rotation again
    [q2new,R2] = Find_Best_Rotation(q1,q2new);
    
    % Compute geodesic distance
    Ec = acos(InnerProd_Q(q1,q2new));
    
    % Update best rotation, re-parameterization, "starting point"
    if Ec < minE            % found better matching of q1 and q2
        Rbest = R2*R;
        q2best = q2new;
        minE = Ec;
        shiftbest = scl*ctr;
        
        if(reparamFlag)
            gamIbest = gamI;
        else
            gamIbest = linspace(0,1,T);
        end
    end
end

return;