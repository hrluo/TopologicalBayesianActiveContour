function c_v = SmoothUpdate(beta,reg)
% Inputs:
% beta = curve
% reg = 1 if want to regularize estimate of curvature (default=0)

% Output:
% curvature = univariate function containing curvature values at each
% point (will be multipled by outward unit normal vector for update)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Settings
N = size(beta,2);
p = 2;      % "pad" curve by two for proper gradient computation

if ~exist('reg','var') || isempty(reg), reg = 0; end
thresh = 0.1;

% Update
[fder,~] = gradient([beta(:,(N-p):(N-1)),beta,beta(:,2:(2+p-1))],1/(N-1+2*p));
[sder,~] = gradient(fder,1/(N-1+2*p));
num = fder(2,:).*sder(1,:)-fder(1,:).*sder(2,:);    % since curve coordinates swap on image
denom = (fder(1,:).^2+fder(2,:).^2).^1.5;
c_v = num./denom;
 
% Regularize against extremely high curvatures by replacing with average
% curvature of two neighboring points
if reg==1
    for i=(p+1):(length(c_v)-p)
        if c_v(i)>thresh
            c_v(i) = thresh;
        end
    end
end

c_v(1:p) = []; c_v((end-p+1):end) = [];
c_v(end) = c_v(1);      % since first point and last point are the same