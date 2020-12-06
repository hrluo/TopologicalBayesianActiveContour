function [JI,JD] = JaccardDist(P,Q)
% Jaccard Distance
% Inputs:
% P = Regions of sets 1. a set of 2-d pixel points representing the whole
% region of the segmented object
% Q = Regions of sets 2. a set of 2-d pixel points representing the whole
% region of the ground truth
 
% Output:
% JI = a numeric value being jaccard index
% JD = a numeric value being jaccard distance (what we wnat to use)

% Ref: https://www.mathworks.com/help/images/ref/jaccard.html if input and
% output are both images.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JD = sum(P & Q)/sum(P | Q);     % Jaccard distance
JI = 1 - JD;                    % Jaccard index
end