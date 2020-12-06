function [seg,totE,interior,tmr] = TOPBACSegT(testI,trainI,trainbeta,n_curves,init,arginit,lambda1,lambda2,lambda3,figs,N_iter,tau,q_bar,bandwidth,sm_init)
%% Perform TOP-BAC segmentation with use of training data
% (e.g., simulations, skin lesion data)

% Inputs:
% testI = image to be segmented (as matrix of pixel values)
% trainI = set of M training images (as matrix of pixel values)
% testI, trainI should have same dimensions
% trainbeta = 2 x N x n_train x n_train_curves matrix representing contours
% from training image - n_train corresponds to number of training images,
% n_train_curves corresponds to number of known contours per training image
% n_curves = number of curves to use for TOP+BAC (may or may not be the
% same as n_train_curves, depending on what user desires)
% init = type of initialization for contour
%   1 = draw contour initializations by hand (default)
%   2 = import initialization mask from external file (must be in current
%   directory) and automatically select n_curves contours with the largest
%   rough estimate of area within
%   3 = import initialization mask from external file (must be in current
%   directory) and cycle through contours until n_curve initializations
%   accepted by user
%   4 = input initialization curve as 2 x N x n_curves matrix (useful for
%   comparing different settings)
% arginit = file name (for init=2,3); 2 x N x n_curves matrix of initial
% contours (for init=6)
% lambda1, lambda2, lambda3 = update constants, either specified as
% n_curve-dim. vector if different constants for each contour desired, or
% as 1 value for each if same constant to be used for all contours, for
% image, smoothing, and shape prior terms, respectively (default 0.3, 0.6,
% 0.05 for all contours)
% figs = 1 if user wants updated plot of contour at each iteration, can
% slow down computation, particularly if image is high resolution
% (default=0)
% N_iter = maximum number of iterations to run the algorithm (default=300)
% tau = cutoff tolerance for algorithm to stop (default=1e-7)
% bandwidth = bandwidth/bin width for pixel density estimator (defaults to
% default of ksdensity/histcounts functions if unspecified)
% sm_init = 1 if user wants to smooth initialized curves, useful if
% initializations are bumpy (default=0)

% Outputs:
% seg = 2 x N x n_curves matrix containing final contours after TOP+BAC
% totE = computed energy at each iteration of algorithm for all n_curves
% interior = binary mask representing interior of seg
% tmr = total time elapsed (in seconds)
%% Basics
[x,y] = size(testI);
[~,N,n_train,n_train_curves] = size(trainbeta);

% Determine if number of curves desired matches number of curves provided
% in training data - dictates estimation of pixel densities, computation of
% prior updates
match = (n_curves==n_train_curves);

% Default parameters
if ~exist('init','var') || isempty(init), init = 1; end
if ~exist('arginit','var') || isempty(arginit), arginit = []; end
if ~exist('lambda1','var') || isempty(lambda1), lambda1 = repmat(0.3,1,n_curves); end
if ~exist('lambda2','var') || isempty(lambda2), lambda2 = repmat(0.6,1,n_curves); end
if ~exist('lambda3','var') || isempty(lambda3), lambda3 = zeros(1,n_curves); end
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
if length(lambda3)==1
    lambda3 = repmat(lambda3,1,n_curves);
end

% If only one initial curve is specified
if (init==4 && (size(arginit,3)~=n_curves))
    tmp = arginit;
    clear arginit
    for i=1:n_curves
        arginit(:,:,i) = tmp;
    end
end

% Indicator of use of shape prior term
for i=1:n_curves
    if lambda3(i)>0
        shapep(i) = 1;
    else
        shapep(i) = 0;
    end
end

%% Compute Karcher mean shape from training curves (if using shape prior)
tic
if any(shapep==1)
    % Compute Karcher mean if not already computed elsewhere
    if ~exist('q_bar','var') || isempty(q_bar)
        for i=1:n_train_curves
            [~,q_bar(:,:,i)] = FindElasticMean(trainbeta(:,:,:,i),1);
        end
    end
    
    % Compute Karcher covariance
    cutoff = 0.8;   % cutoff for percent total variation for PCA
    for i=1:n_train_curves
        for j=1:n_train
            q_train(:,:,j,i) = curve_to_q(trainbeta(:,:,j,i));
        end
        K = FindElasticCovariance(q_bar(:,:,i),q_train(:,:,:,i));
        [U,S,~] = svd(K);
        CV = cumsum(diag(S(1:(n_train-1),1:(n_train-1))))/sum(diag(S(1:(n_train-1),1:(n_train-1))));
        m(i) = find(CV>=cutoff,1);
        Um{i} = U(:,1:m(i));
        Sm{i} = S(1:m(i),1:m(i));
        Sm_inv{i} = zeros(m(i),m(i));
        for j=1:m(i)
            Sm_inv{i}(j,j) = 1/Sm{i}(j,j);
        end
        
        % Used for prior gradient update
        delta(i) = 0.8*(1/Sm_inv{i}(m(i),m(i)));    % chosen to be smaller than smallest kept eigenvalue
        A{i} = Um{i}*Sm_inv{i}*Um{i}' + (1/delta(i)^2)*(eye(2*N,2*N)-Um{i}*Um{i}');
    end
end

%% Initialize contour
if init==1  % Draw contour by hand
    % Query user to draw polygon initialization - click to introduce a new
    % point, and then double click once at the end to finish
    % NOTE: Specify polygon in a clockwise orientation to ensure outward
    % unit normal is defined correctly, as well as to match orientation of
    % prior.
    for i=1:n_curves
        [~,tmp(2,:),tmp(1,:)] = roipoly(testI);
        arginit(:,:,i) = ReSampleCurve(tmp,N);
        clear tmp
    end
    
elseif init==2  % Import TOP result (as mask) and automatically select first n_curve largest boundary
    % curves by area as initializations.
    cd('TOP Initializations')
    seg_mask = imread(arginit);
    clear arginit
    cd('../')
    [boundary,~,K] = bwboundaries(seg_mask(:,:,1)); % K = number of objects detected in image
    
    % Truncate boundary cell to first K objects and compute area contained
    % within each boundary, sorting from largest to smallest
    new_boundary = boundary(1:K);
    for i=1:K
        ArI(i) = polyarea(new_boundary{i}(:,2)',new_boundary{i}(:,1)');
    end
    
    [~,idx] = sort(ArI,'descend');
    
    % Resample contour to N points (to match length of training curve)
    for i=1:n_curves
        arginit(:,:,i) = ReSampleCurve(new_boundary{idx(i)}',N);
    end
    
elseif init==3 % Import TOP result (as mask) and cycle through boundary curves until n_curves are selected.
    cd('TOP Initializations')
    seg_mask = imread(arginit);
    clear arginit
    cd('../')
    [boundary,~,K] = bwboundaries(seg_mask(:,:,1)); % K = number of objects detected in image
    
    % Truncate boundary cell to first K objects and compute area contained
    % within each boundary, sorting from largest to smallest
    new_boundary = boundary(1:K);
    for i=1:K
        ArI(i) = polyarea(new_boundary{i}(:,2)',new_boundary{i}(:,1)');
    end
    
    [~,idx] = sort(ArI,'descend');
    
    % Cycle through boundary curves until n_curves accepted
    n_acc = 0;
    
    for i=1:K
        tmp = new_boundary{idx(i)}';
        
        % Plot ith boundary curve
        imshow(testI)
        hold on
        plot(tmp(2,:),tmp(1,:),'r','LineWidth',3)
        
        % Query user to accept or reject the boundary curve
        prompt = 'Would you like to proceed with this initialization? (1 = yes, 0 = no) ';
        u_tmp = input(prompt);
        if u_tmp==1
            arginit(:,:,n_acc+1) = ReSampleCurve(tmp,N);
            n_acc = n_acc+1;
        end
        
        clf
        
        % Stop when n_acc = n_curves+1;
        if n_acc==n_curves
            break
        end
    end
    
    % If user exhausts list without selecting enough curves, re-define
    % n_curves to be the total number accepted
    if n_acc < n_curves
        n_curves = n_acc;
    end
    
    close all
    
else           % Use input initialization curve
    tmp = arginit;
    clear arginit
    
    % Resample contour to N points (to match length of training curve)
    for i=1:n_curves
        arginit(:,:,i) = ReSampleCurve(tmp(:,:,i),N);
    end
end

% Check to see if the initialization is appropriate for the image
if figs==1
    imshow(testI)
    hold on
    for i=1:n_curves
        plot(arginit(2,:,i),arginit(1,:,i),'r','LineWidth',3)
    end
end

% Pre-smooth TOP initialized contours (since some of these are extremely
% rough)
if sm_init==1
    span = 11;  % span for moving average smoother
    w = (span-1)/2;
    for i=1:n_curves
        % Closed contours - need to pad with neighboring points for smoother
        tmp = [arginit(:,(end-w):(end-1),i),arginit(:,:,i),arginit(:,2:(1+w),i)];
        initc{i}(1,:) = smooth(tmp(1,:),span);
        initc{i}(2,:) = smooth(tmp(2,:),span);
        initc{i}(:,[1:w,(end-w+1):end]) = [];
        initc{i}(:,end) = initc{i}(:,1);
    end
else
    for i=1:n_curves
        initc{i} = arginit(:,:,i);
    end
end

%% Estimate interior/exterior pixel intensity densities from training images/curves
M = 255;    % number of pixel values between [0,1] to define densities at

% Convert images to double format if not already in it
if ~isa(testI,'double')
    testI = im2double(testI);
end
testI = testI(:,:,1);  % in case image is greyscale

for i=1:n_train
    if ~isa(trainI{i},'double')
        trainI{i} = im2double(trainI{i});
    end
    trainI{i} = trainI{i}(:,:,1);  % in case image is greyscale
end

% Density estimators of pixel values interior and exterior to ground truth
for i=1:n_train_curves
    [p_in{i},p_out{i},llik{i}] = TrainingPixelDensity(trainI,trainbeta(:,:,:,i),M,bandwidth);
    nlp_in{i} = -log(p_in{i});    % negative log likelihood for interior pixel values
    nlp_out{i} = -log(p_out{i});  % negative log likelihood for exterior pixel values
    
    % Compute log-likelihood on each pixel value of test image
    for j=1:x
        for k=1:y
            pix = testI(j,k);
            [~,idx] = min(abs(repmat(pix,1,M)-linspace(0,1,M)),[],2);
            llikI{i}(j,k) = llik{i}(idx);
        end
    end
end

%% Compute initial energy values
for i=1:n_curves
    % If same number of curves desired as input for training data, then use
    % estimates from each curve separately; otherwise, use estimates from
    % the 1 contour.
    if match
        idx = i;
    else
        idx = 1;
    end
    
    % Image energy term
    E_img{i}(1) = ImageEnergy(testI,initc{i},nlp_in{idx},nlp_out{idx});
    
    % Smoothing energy term
    E_smooth{i}(1) = SmoothEnergy(initc{i});
    
    % Total energy
    totE{i}(1) = E_img{i}(1)+E_smooth{i}(1);
    
    % Prior energy term (and sum to total)
    if shapep(i)==1
        [E_prior{i}(1),d_prior{i}(1)] = PriorEnergy(initc{i},q_bar(:,:,idx),Um{idx},Sm_inv{idx},delta(idx));
    else
        E_prior{i}(1) = 0;
        d_prior{i}(1) = 0;
    end
    totE{i}(1) = totE{i}(1)+E_prior{i}(1);
end

%% Update segmentation and iterate (separate BAC for each curve)
for j=1:n_curves
    if match
        idx = j;
    else
        idx = 1;
    end
    
    % Only call from cells once to speed up computation
    tmp_init = initc{j};
    tmp = tmp_init;
    tmp_totE = totE{j};
    tmp_Eimg = E_img{j};
    tmp_Esmooth = E_smooth{j};
    tmp_Eprior = E_prior{j};
    tmp_llikI = llikI{idx};
    tmp_nlp_in = nlp_in{idx};
    tmp_nlp_out = nlp_out{idx};
    if shapep(j)==1
        tmp_qbar = q_bar(:,:,idx);
        tmp_A = A{idx};
        tmp_Um = Um{idx};
        tmp_Sminv = Sm_inv{idx};
        tmp_delta = delta(idx);
    end
    
    for i=2:N_iter
        %% Prior update
        if shapep(j)==1
            [delP,tmp] = PriorUpdate(tmp,tmp_qbar,tmp_A,0);
        end
        
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
        if shapep(j)==1
            if lambda1(j)==0 && lambda2(j)>0
                tmp = tmp-lambda2(j)*delS-lambda3(j)*delP;
            elseif lambda2(j)==0 && lambda1(j)>0
                tmp = tmp-lambda1(j)*delI-lambda3(j)*delP;
            elseif lambda1(j)==0 && lambda2(j)==0
                tmp = tmp-lambda3(j)*delP;
            else
                tmp = tmp-lambda1(j)*delI-lambda2(j)*delS-lambda3(j)*delP;
            end
        else
            if lambda1(j)==0
                tmp = tmp-lambda2(j)*delS;
            elseif lambda2(j)==0
                tmp = tmp-lambda1(j)*delI;
            else
                tmp = tmp-lambda1(j)*delI-lambda2(j)*delS;
            end
        end
        tmp(:,end) = tmp(:,1);
        tmp = ReSampleCurve(tmp,N);
        
        %% Updated energy values
        tmp_Eimg(i) = ImageEnergy(testI,tmp,tmp_nlp_in,tmp_nlp_out);
        tmp_Esmooth(i) = SmoothEnergy(tmp);
        tmp_totE(i) = tmp_Eimg(i)+tmp_Esmooth(i);
        if shapep(j)==1
            [tmp_Eprior(i),tmp_dprior(i)] = PriorEnergy(tmp,tmp_qbar,tmp_Um,tmp_Sminv,tmp_delta);
        else
            tmp_Eprior(i) = 0;
            tmp_dprior(i) = 0;
        end
        tmp_totE(i) = tmp_totE(i)+tmp_Eprior(i);
        tmp_smooth = smooth(tmp_totE);
        sc = sum(diff(tmp_smooth(max(1,i-5):i))>0);
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
    E_prior{j} = tmp_Eprior;
    d_prior{j} = tmp_dprior;
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