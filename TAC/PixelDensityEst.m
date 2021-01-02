function [p_in,p_out,llik] = PixelDensityEst(I,beta,s,bandwidth)
% Inputs:
% I = test image (as cell) - for neuron example, estimate interior and
% exterior pixel density based on the TOP initialization as a crude
% estimator
% beta = 2 x N x n_curves vector-valued curve representing TOP
% initializations for test image
% s = number of equally-spaced pixel values between 0 and 1
% bandwidth = bandwidth for kernel density estimator or histogram estimator

% Outputs:
% p_in = estimated density of interior pixel values
% p_out = estimated density of exterior pixel values
% llik = log(p_in/p_out) = log likelihood for pixel value memberships
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,~,n_curves] = size(beta);

% Augment all interior and exterior pixel values from image using all of
% the TOP initialization curves
for i=1:n_curves
    pmask(:,:,i) = roipoly(I,beta(2,:,i),beta(1,:,i));
end

mask = (sum(pmask,3)>0);  % "cumulative" mask
vmaski = mask(:);
vimg = I(:);
vmasko = logical(abs(vmaski-1));
pix_in = vimg(vmaski);
pix_out = vimg(vmasko);

% Estimate densities for pixel value space
if any(pix_in > 0 & pix_in < 1) || any(pix_out > 0 & pix_out < 1)
    % Kernel density estimator if pixel values are defined between 0 and 1
    eps = 1e-10;    % to ensure density is fully specified over [0,1]
    if ~exist('bandwidth','var') || isempty(bandwidth)
        [p_in,~,bandwidth] = ksdensity(pix_in,linspace(0,1,s),'BoundaryCorrection','reflection','Support',[-eps 1+eps]);
        p_out = ksdensity(pix_out,linspace(0,1,s),'Bandwidth',bandwidth,'BoundaryCorrection','reflection','Support',[-eps 1+eps]);
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