function [p_in,p_out,llik] = TrainingPixelDensity(I,beta,s,bandwidth)
% Inputs:
% I = set of n_train training images (as a cell array, each entry is a
% separate image) - implicitly, we assume that these images all have
% similar distribution of pixel densities
% beta = 2 x N x n_train vector-valued curve representing contour from
% training image
% s = number of equally-spaced pixel values between 0 and 1
% bandwidth = bandwidth for kernel density estimator or histogram estimator

% Outputs:
% p_in = estimated density of pixel values inside training image contours
% p_out = estimated density of pixel values outside training image contours
% llik = log(p_in/p_out) = log likelihood for pixel value memberships
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_train = length(I);    % number of training images

% Augment all interior and exterior pixel values from all training images
pix_in = [];
pix_out = [];

for i=1:n_train
    mask = roipoly(I{i},beta(2,:,i),beta(1,:,i));
    vmaski = mask(:);
    vimg = I{i}(:);
    vmasko = logical(abs(vmaski-1));
    pix_in = [pix_in;vimg(vmaski)];
    pix_out = [pix_out;vimg(vmasko)];
end

% Estimate densities for pixel value space
if any(pix_in > 0 & pix_in < 1) || any(pix_out > 0 & pix_out < 1)
    % Kernel density estimator if pixel values are defined between 0 and 1
    eps = 1e-10;    % to ensure density is fully specified over [0,1]
    if ~exist('bandwidth','var') || isempty(bandwidth)
        p_in = ksdensity(pix_in,linspace(0,1,s),'BoundaryCorrection','reflection','Support',[-eps 1+eps]);
        p_out = ksdensity(pix_out,linspace(0,1,s),'BoundaryCorrection','reflection','Support',[-eps 1+eps]);
    else
        p_in = ksdensity(pix_in,linspace(0,1,s),'Bandwidth',bandwidth,'BoundaryCorrection','reflection','Support',[-eps 1+eps]);
        p_out = ksdensity(pix_out,linspace(0,1,s),'Bandwidth',bandwidth,'BoundaryCorrection','reflection','Support',[-eps 1+eps]);
    end
else
    % Histogram estimator if pixel values are binary (either 0 or 1)
    if ~exist('bandwidth','var') || isempty(bandwidth)
        [p_in,~] = histcounts(pix_in,255,'BinLimits',[0 1],'Normalization','pdf');
        [p_out,~] = histcounts(pix_out,255,'BinLimits',[0 1],'Normalization','pdf');
    else
        [p_in,~] = histcounts(pix_in,255,'BinWidth',bandwidth,'BinLimits',[0 1],'Normalization','pdf');
        [p_out,~] = histcounts(pix_out,255,'BinWidth',bandwidth,'BinLimits',[0 1],'Normalization','pdf');
    end
end

% Compute log-likelihood
llik = log(p_in./p_out);
llik(isnan(llik)) = 0;
for i=1:s
    if isinf(llik(i))
        if llik(i)<0
            p_in(i) = 1e-100;
            llik(i) = log(p_in(i)/p_out(i));
        else
            p_out(i) = 1e-100;
            llik(i) = log(p_in(i)/p_out(i));
        end
    end
end