function [HPQ,HQP,HD] = HammingDist(P,Q)
% Hamming Distance
% Inputs:
% P = Regions of set 1 - a set of 2-d pixel points representing the whole
% region of the segmented object
% Q = Regions of set 2 - a set of 2-d pixel points representing the whole
% region of the ground truth

% For definition, see "Performance Evaluation of Image Segmentation" by 
% Fernando C. Monteiro and Aurlio C. Campilho

% Output:
% HPQ = Hamming P->Q distance
% HQP = Hamming Q->P distance
% HD = Hamming distance (what we want to use)

% This distance assumes that we have only one object in the image,
% represented by values 1.
% In image P, there are only two clusters: those pixel points with value 1;
% those pixel points with value 0; Similar in image Q.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sizevec = size(P);
img_size = sizevec(1)*sizevec(2);
cluster_P_1 = P;
cluster_P_2 = 1-P;
cluster_Q_1 = Q;
cluster_Q_2 = 1-Q;

% N = nnz(X) returns the number of nonzero elements in matrix X
HPQ = nnz(cluster_P_1 & cluster_Q_2)+nnz(cluster_P_2 & cluster_Q_1);
HQP = nnz(cluster_Q_1 & cluster_P_2)+nnz(cluster_Q_2 & cluster_P_1);
HD = (HPQ+HQP)/(2*img_size);
end