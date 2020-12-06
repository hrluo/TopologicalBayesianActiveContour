function [seg,totE,interior,tmr] = TOPBACSeg(testI,trainI,trainbeta,n_curves,init,arginit,lambda1,lambda2,lambda3,figs,N_iter,tau,q_bar,bandwidth)
% Inputs:
% testI = image to be segmented (as matrix of pixel values)
% trainI = set of M training images (as matrix of pixel values)
% testI, trainI should have same dimensions
% trainbeta = 2 x N x M matrix representing contours from training image
% (if only a single contour available from all training data)
% n_curves = number of curves to use for TOP+BAC
% init = type of initialization for contour
%   1 = draw contour initializations by hand (default)
%   2 = import initialization mask from external file (must be in current
%   directory) and automatically select n_curve contours with the largest
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

%% Estimate interior/exterior pixel intensity densities from training images/curves
M = 255;    % number of pixel values between [0,1] to define densities at

% Convert images to double format if not already in it
if ~isa(testI,'double')
    testI = im2double(testI);
end
for i=1:n_train
    if ~isa(trainI{i},'double')
        trainI{i} = im2double(trainI{i});
    end
end

% Density estimators of pixel values interior and exterior to ground truth
for i=1:n_train_curves
    [p_in{i},p_out{i},llik{i}] = TrainingPixelDensity(trainI,trainbeta(:,:,:,i),M,bandwidth);
    nlp_in{i} = -log(p_in{i});    % negative log likelihood for interior pixel values
    nlp_out{i} = -log(p_out{i});  % negative log likelihood for exterior pixel values
end

% Compute log-likelihood on each pixel value of test image
for p=1:n_train_curves
    for i=1:x
        for j=1:y
            pix = testI(i,j);
            [~,idx] = min(abs(repmat(pix,1,M)-linspace(0,1,M)),[],2);
            llikI{p}(i,j) = llik{p}(idx);
        end
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
        seg(:,:,i) = ReSampleCurve(tmp,N);
        clear tmp
    end
    
elseif init==2  % Import TOP result (as mask) and automatically select first n_curve largest boundary
                % curves by area as initializations.
    cd('TOP Initializations')
    seg_mask = imread(arginit);
    cd('../')
    [boundary,~,K] = bwboundaries(seg_mask(:,:,1)); % K = number of objects detected in image
    
    % Truncate boundary cell to first K objects and compute rough estimate
    % of area contained within each boundary
    new_boundary = boundary(1:K);
    for i=1:K
        rng = range(new_boundary{i});
        Ar(i) = rng(1)*rng(2);
    end
    
    % Select boundary curve that has the highest number of sampling points
    [~,idx] = sort(Ar,'descend');

    % Resample contour to N points (to match length of training curve)
    for i=1:n_curves
        seg(:,:,i) = ReSampleCurve(new_boundary{idx(i)}',N);
    end
    
elseif init==3 % Import TOP result (as mask) and cycle through boundary curves until n_curves are selected.
    cd('TOP')
    seg_mask = imread(arginit);
    cd('../')
    [boundary,~,K] = bwboundaries(seg_mask(:,:,1)); % K = number of objects detected in image
    
    % Truncate boundary cell to first K objects and compute rough estimate
    % of area contained within each boundary
    new_boundary = boundary(1:K);
    for i=1:K
        rng = range(new_boundary{i});
        Ar(i) = rng(1)*rng(2);
    end
    
    % Select boundary curve that has the highest number of sampling points
    [~,idx] = sort(Ar,'descend');

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
            seg(:,:,n_acc+1) = ReSampleCurve(tmp,N);
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
    % Resample contour to N points (to match length of training curve)
    for i=1:n_curves
        seg(:,:,i) = ReSampleCurve(arginit(:,:,i),N);
    end
end

% Check to see if the initialization is appropriate for the image
if figs==1
    imshow(testI)
    hold on
    for i=1:n_curves
        plot(seg(2,:,i),seg(1,:,i),'r','LineWidth',3)
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
    
    E_img(1,i) = ImageEnergy(testI,seg(:,:,i),nlp_in{idx},nlp_out{idx});
    E_smooth(1,i) = SmoothEnergy(seg(:,:,i));
    totE(1,i) = E_img(1,i)+E_smooth(1,i);

    if shapep(i)==1
        [E_prior(1,i),d_prior(1,i)] = PriorEnergy(seg(:,:,i),q_bar(:,:,idx),Um{idx},Sm_inv{idx},delta(idx));
        totE(1,i) = totE(1,i)+E_prior(1,i);
    end
end

% Indicator of whether or not a particular contour is active (i.e., still
% being updated or has reached convergence) - active if = 0
ind = zeros(1,n_curves);

% Initialize smoothtotE
smoothtotE = zeros(N_iter,n_curves);

%% Update segmentation and iterate
for i=2:N_iter
    for j=1:n_curves
        % If indicator is 0, proceed with update for jth curve; otherwise,
        % move to next curve in loop
        if ind(j)==0
            if match
                idx = j;
            else
                idx = 1;
            end
            
            %% Shape prior gradient update
            if shapep(j)==1
                [delP,seg(:,:,j)] = PriorUpdate(seg(:,:,j),q_bar(:,:,idx),A{idx},0);
            end
            
            %% Compute outward unit normal vector
            nrm_v = OutwardUnitNormal(seg(:,:,j));
            
            %% Image gradient update
            sc_v = ImageUpdate(seg(:,:,j),llikI{idx});
            delI(1,:) = sc_v.*nrm_v(1,:);
            delI(2,:) = sc_v.*nrm_v(2,:);
            
            %% Smoothing gradient update
            c_v = SmoothUpdate(seg(:,:,j));
            delS(1,:) = c_v.*nrm_v(1,:);
            delS(2,:) = c_v.*nrm_v(2,:);
            
            %% Combine all for full update
            if shapep(j)==1
                if lambda1(j)==0 && lambda2(j)>0
                    seg(:,:,j) = seg(:,:,j)-lambda2(j)*delS-lambda3(j)*delP;
                elseif lambda2(j)==0 && lambda1(j)>0
                    seg(:,:,j) = seg(:,:,j)-lambda1(j)*delI-lambda3(j)*delP;
                elseif lambda1(j)==0 && lambda2(j)==0
                    seg(:,:,j) = seg(:,:,j)-lambda3(j)*delP;
                else
                    seg(:,:,j) = seg(:,:,j)-lambda1(j)*delI-lambda2(j)*delS-lambda3(j)*delP;
                end
            else
                if lambda1(j)==0
                    seg(:,:,j) = seg(:,:,j)-lambda2(j)*delS;
                elseif lambda2(j)==0
                    seg(:,:,j) = seg(:,:,j)-lambda1(j)*delI;
                else
                    seg(:,:,j) = seg(:,:,j)-lambda1(j)*delI-lambda2(j)*delS;
                end
            end
            seg(:,end,j) = seg(:,1,j);
            seg(:,:,j) = ReSampleCurve(seg(:,:,j),N);
        end
    end
    
    if figs==1
        figure(1)
        imshow(testI)
        hold on
        for j=1:n_curves
            plot(seg(2,:,j),seg(1,:,j),'r','LineWidth',3)
        end
    end
    
    %% Update energy values
    for j=1:n_curves
        if ind(j)==1        % inactive contour
            E_img(i,j) = E_img(i-1,j);
            E_smooth(i,j) = E_smooth(i-1,j);
            if shapep(j)==1
                E_prior(i,j) = E_prior(i-1,j);
                d_prior(i,j) = d_prior(i-1,j);
            end
            totE(i,j) = totE(i-1,j);
            smoothtotE(i,j) = smoothtotE(i-1,j);
        else
            if match
                idx = j;
            else
                idx = 1;
            end
            
            E_img(i,j) = ImageEnergy(testI,seg(:,:,j),nlp_in{idx},nlp_out{idx});
            E_smooth(i,j) = SmoothEnergy(seg(:,:,j));
            totE(i,j) = E_img(i,j)+E_smooth(i,j);
            if shapep(j)==1
                [E_prior(i,j),d_prior(i,j)] = PriorEnergy(seg(:,:,j),q_bar(:,:,idx),Um{idx},Sm_inv{idx},delta(idx));
                totE(i,j) = totE(i,j)+E_prior(i,j);
            end
            
            % Smooth total energy for stopping criterion
            smoothtotE(1:i,j) = smooth(totE(:,j));
            
            % Stop algorithm for jth curve if its total energy converges or 3
            % consecutive steps result in increased energy
            if (abs((smoothtotE(i,j)-smoothtotE(i-1,j)))/abs(smoothtotE(i,j)) < tau) || (i>4 && (smoothtotE(i,j)>smoothtotE(i-1,j)) && (smoothtotE(i-1,j)>smoothtotE(i-2,j)) && (smoothtotE(i-2,j)>smoothtotE(i-3,j)))
                ind(j) = 1;
            else
                ind(j) = 0;
            end
        end
    end
    
    % Exit loop if all contours have converged
    if sum(ind)==n_curves
        break
    end
    
    i
end

tmr = toc;

% Binary mask indicating interior of final contours
interior = zeros(x,y);
for i=1:n_curves
    interior_c(:,:,i) = roipoly(testI,seg(2,:,i),seg(1,:,i));
    interior = interior+interior_c(:,:,i);
end
interior = mod(interior,2);