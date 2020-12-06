function [HauD,JD,HamD,PM,ESD] = SegDistTop(estInt,trueInt,estBeta,trueBeta)
% Inputs:
% estInt = estimated segmentation as a binary image mask
% trueInt = ground truth segmentation as a binary image mask
% estBeta = estimated contour of boundary as 2 x N x C matrix, where C is
% the number of curves
% trueBeta = ground truth contour of boundary as 2 x N x C matrix, where C
% is the number of curves

% Outputs:
% HauD = Hausdorff distance
% JD = Jaccard distance
% HamD = Hamming distance
% PM = performance measure distance
% ESD = re-scaled elastic shape distance for each curve (bounded by pi/2?)

% Note: This can be used for images with just 1 connected component or more
% than 1 connected component.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[HauD,~]= HausdorffDist(estInt,trueInt,1);
[~,JD] = JaccardDist(estInt,trueInt);
[~,~,HamD] = HammingDist(estInt,trueInt);

% Performance measure distance
TP = sum(sum(estInt & trueInt));         % True positive
TN = sum(sum((1-estInt) & (1-trueInt))); % True negative
FP = sum(sum(estInt & (1-trueInt)));     % False positive
FN = sum(sum((1-estInt) & trueInt));     % False negative
PM = 1-(TP)/(TP+FN+FP);

% Elastic shape distance for each curve
n_curves = size(estBeta,3);
for i=1:n_curves
    est_q = curve_to_q(estBeta(:,:,i));
    true_q = curve_to_q(trueBeta(:,:,i));
    [true_qn,~,~,~] = Find_Rotation_and_Seed_unique(est_q,true_q,1);
    ESD(i) = 2*acos(InnerProd_Q(est_q,true_qn))/pi;
end