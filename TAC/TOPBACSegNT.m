function [seg,totE,interior,tmr] = TOPBACSegNT(testI,arginit,lambda1,lambda2,figs,N_iter,tau,bandwidth,sm_init,pool)
%% Perform TOP-BAC segmentation without use of training data
% (e.g., neuron data)

% Inputs:
% testI = image to be segmented (as matrix of pixel values)
% arginit = 2 x N x n_curves dimensional initial contour
% lambda1, lambda2 = update constants either specified as
% n_curve-dim. vector if different constants for each contour desired, or
% as 1 value for each if same constant to be used for all contours, for
% image and smoothing terms, respectively (default 0.3, 0.6 for all
% contours)
% figs = 1 if user wants updated plot of contour at each iteration, can
% slow down computation, particularly if image is high resolution
% (default=0)
% N_iter = maximum number of iterations to run the algorithm (default=300)
% tau = cutoff tolerance for algorithm to stop (default=1e-7)
% bandwidth = bandwidth/bin width for pixel density estimator (defaults to
% default of ksdensity/histcounts functions if unspecified)
% sm_init = 1 if user wants to smooth initialized curves, useful if
% initializations are bumpy (default=0)
% pool = 1 if want to pool all contours together to estimate interior and
% exterior pixel densities (default=0)

% Output:
% seg = 2 x N x n_curves matrix containing final contours after TOP+BAC
% totE = computed energy at each iteration of algorithm for all n_curves
% interior = binary mask representing interior of seg
% tmr = total time elapsed (in seconds)
%% Basics
tic
[x,y,~] = size(testI);
[~,N,n_curves] = size(arginit);

% Default parameters
if ~exist('lambda1','var') || isempty(lambda1), lambda1 = 0.3; end
if ~exist('lambda2','var') || isempty(lambda2), lambda2 = 0.6; end
if ~exist('figs','var') || isempty(figs), figs = 0; end
if ~exist('N_iter','var') || isempty(N_iter), N_iter = 300; end
if ~exist('tau','var') || isempty(tau), tau = 1e-7; end
if ~exist('bandwidth','var') || isempty(bandwidth), bandwidth = []; end
if ~exist('sm_init','var') || isempty(sm_init), sm_init = 0; end
if ~exist('pool','var') || isempty(pool), pool = 0; end

% If only one value specified for lambda1, lambda2, lambda3, then create
% vectors of size n_curves with repeated values
if length(lambda1)==1
    lambda1 = repmat(lambda1,1,n_curves);
end
if length(lambda2)==1
    lambda2 = repmat(lambda2,1,n_curves);
end

% Pre-smooth TOP initialized contours (since some of these are extremely
% rough)
if sm_init==1
    span = 11;  % span for moving average smoother
    w = (span-1)/2;
    for i=1:n_curves
        % Closed contours - need to pad with neighboring points for smoother
        tmp = [arginit(:,(end-w):(end-1),i),arginit(:,:,i),arginit(:,2:(1+w),i)];
        init{i}(1,:) = smooth(tmp(1,:),span);
        init{i}(2,:) = smooth(tmp(2,:),span);
        init{i}(:,[1:w,(end-w+1):end]) = [];
        init{i}(:,end) = init{i}(:,1);
    end
else
    for i=1:n_curves
        init{i} = arginit(:,:,i);
    end
end

%% Estimate interior/exterior pixel intensity densities from training images/curves
M = 255;    % number of pixel values between [0,1] to define densities at

% Convert image to double format if not already in it
if ~isa(testI,'double')
    testI = im2double(testI);
end
testI = testI(:,:,1);  % in case image is greyscale

% Pixel densities of interior and exterior estimated from TOP contour
% initializations (since we do not have training data)
if pool==1  % estimate one set of densities based on union of initializations
    [tmp_in,tmp_out,tmp_llik] = PixelDensityEst(testI,arginit,M,bandwidth);
    tmp_nlp_in = -log(tmp_in);
    tmp_nlp_out = -log(tmp_out);
    
    % Compute log-likelihood on each pixel value of test image
    for j=1:x
        for k=1:y
            pix = testI(j,k);
            [~,idx] = min(abs(repmat(pix,1,M)-linspace(0,1,M)),[],2);
            tmp_llikI(j,k) = tmp_llik(idx);
        end
    end
    
    % Save a copy for each curve
    for i=1:n_curves
        p_in{i} = tmp_in;
        p_out{i} = tmp_out;
        nlp_in{i} = tmp_nlp_in;
        nlp_out{i} = tmp_nlp_out;
        llik{i} = tmp_llik;
        llikI{i} = tmp_llikI;
    end
else  % estimate separate densities for each initialization
    for i=1:n_curves
        I{1} = testI;
        [p_in{i},p_out{i},llik{i}] = TrainingPixelDensity(I,arginit(:,:,i),M,bandwidth);
        nlp_in{i} = -log(p_in{i});
        nlp_out{i} = -log(p_out{i});
        
        % Compute log-likelihood on each pixel value of test image
        for j=1:x
            for k=1:y
                pix = testI(j,k);
                [~,idx] = min(abs(repmat(pix,1,M)-linspace(0,1,M)),[],2);
                llikI{i}(j,k) = llik{i}(idx);
            end
        end
        
        i
    end
end

%% Compute initial energy values
for i=1:n_curves
    % Image energy term
    E_img{i}(1) = ImageEnergy(testI,init{i},nlp_in{i},nlp_out{i});

    % Smoothing energy term
    E_smooth{i}(1) = SmoothEnergy(init{i});
    
    % Total energy
    totE{i}(1) = E_img{i}(1)+E_smooth{i}(1);
end

%% Update segmentation and iterate (separate BAC for each curve)
for j=35:-1:1
    % Only call from cells once to speed up computation
    tmp_init = init{j};
    tmp = tmp_init;
    tmp_totE = totE{j};
    tmp_Eimg = E_img{j};
    tmp_Esmooth = E_smooth{j};
    tmp_llikI = llikI{j};
    tmp_nlp_in = nlp_in{j};
    tmp_nlp_out = nlp_out{j};
    
    for i=2:N_iter
        %% Compute outward unit normal vector
        nrm_v = OutwardUnitNormal(tmp);
        
        %% Image gradient update
        sc_v = ImageUpdate(tmp,tmp_llikI);
        delI(1,:) = sc_v.*nrm_v(1,:);
        delI(2,:) = sc_v.*nrm_v(2,:);
        
        %% Smoothing gradient update
        c_v = SmoothUpdate(tmp,1);
        delS(1,:) = c_v.*nrm_v(1,:);
        delS(2,:) = c_v.*nrm_v(2,:);
        
        %% Full update
        if lambda1(j)==0
            tmp = tmp-lambda2(j)*delS;
        elseif lambda2(j)==0
            tmp = tmp-lambda1(j)*delI;
        else
            tmp = tmp-lambda1(j)*delI-lambda2(j)*delS;
        end
        tmp(:,end) = tmp(:,1);
        tmp = ReSampleCurve(tmp,N);
        
        %% Updated energy values
        tmp_Eimg(i) = ImageEnergy(testI,tmp,tmp_nlp_in,tmp_nlp_out);
        tmp_Esmooth(i) = SmoothEnergy(tmp);
        tmp_totE(i) = tmp_Eimg(i)+tmp_Esmooth(i);
        tmp_smooth = smooth(tmp_totE);
        sc = sum(diff(tmp_smooth(max(1,i-4):i))>0);
        % used to stop algorithm if smoothed energy increasing for 5
        % successive iterations
        
        if (figs==1)
            imshow(testI)
            hold on
            plot(tmp_init(2,:),tmp_init(1,:),'b','LineWidth',3)
            plot(tmp(2,:),tmp(1,:),'r','LineWidth',3)
        end
        
        %% Stopping criterion
        if (abs((tmp_smooth(i)-tmp_smooth(i-1)))/abs(tmp_smooth(i)) < tau) || sc==5
            break
        end
        
        [j,i]
    end
    
    % Save to cells
    seg2{j} = tmp;
    totE{j} = tmp_totE;
    E_img{j} = tmp_Eimg;
    E_smooth{j} = tmp_Esmooth;
end

tmr = toc;

%% Display curves
if figs==1
    figure(1)
    imshow(testI)
    hold on
    for j=1:n_curves
        plot(init{j}(2,:),init{j}(1,:),'b','LineWidth',3)
        plot(seg2{j}(2,:),seg2{j}(1,:),'r','LineWidth',3)
    end
end

%% Binary mask indicating interior of final contours
for i=1:n_curves
    seg(:,:,i) = seg2{i};
end

interior = zeros(x,y);
for i=1:n_curves
    interior_c(:,:,i) = roipoly(testI,seg(2,:,i),seg(1,:,i));
    interior = interior+interior_c(:,:,i);
end
interior = mod(interior,2);
interior = uint8(interior);