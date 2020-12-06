%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for neuron data example %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths - start in D_home folder
D_home = '/Users/jdstrait/OneDrive - University of Georgia/Research/TOPBAC/COMPARE/REVISION/updated_code';
D_data = 'Data/Neuron';

%% Settings
N = 200;  % number of points to resample curves to

%% Extract TOP initialization
% Load in image and full representation of TOP map
cd(D_data)
I = imread('0085.ppm');  % image to be segmented
TOPI = imread('0085_full_2_4_1_add5.ppm');  % TOP segmentation
cd(D_home)

% Find unique RGB pixel values
TOPI_rs = reshape(TOPI,[],3);
upix = unique(TOPI_rs,'rows');

% Form binary mask for each unique RGB pixel value, and extract largest
% boundary from each
for i=1:size(upix,1)
    % Binary mask
    bin_mask{i} = (TOPI(:,:,1)==upix(i,1)) & (TOPI(:,:,2)==upix(i,2)) & (TOPI(:,:,3)==upix(i,3));

    % Extract largest boundary
    clear tmp K new_boundary ArI_can idx
    [tmp,~,K] = bwboundaries(bin_mask{i});
    new_boundary = tmp(1:K);
    for j=1:K
        ArI_can(j) = polyarea(new_boundary{j}(:,2)',new_boundary{j}(:,1)');
    end
    [~,idx] = sort(ArI_can,'descend');
    ArI(i) = ArI_can(idx(1));
    tmp_boundary{i} = new_boundary{idx(1)};
    
    % Compute its centroid, and distance of centroid from the middle pixel
    dc(i) = norm((mean(tmp_boundary{i})-0.5*size(TOPI,1:2))./(size(TOPI,1:2)));
end

% Narrow down to closest to center of image
[~,idx2] = sort(dc,'ascend');
tmp2_boundary = tmp_boundary(idx2(1:ceil(0.50*length(idx2))));
Ar2I = ArI(idx2(1:ceil(0.50*length(idx2))));

% Remove the smallest neurons
[~,idx3] = sort(Ar2I,'descend');
tmp3_boundary = tmp2_boundary(idx3(1:ceil(0.40*length(idx3))));
Ar3I = Ar2I(idx3(1:ceil(0.40*length(idx3))));

% Remove any outliers with respect to size after removing smallest
tmp4_boundary = tmp3_boundary(~isoutlier(Ar3I,'mean'));
Ar4I = Ar3I(~isoutlier(Ar3I,'mean'));

% Compare each boundary to a circle via elastic shape distance to select
% the ones that are "most circular"
theta = linspace(1,0,100);
circ = [cos(2*pi*theta);sin(2*pi*theta)];
q_circ = curve_to_q(circ);

for i=1:length(Ar4I)
    init_TOPI(:,:,i) = ReSampleCurve(tmp4_boundary{i}',N);
    tmp_beta = ReSampleCurve(init_TOPI(:,:,i),100);
    q_tmp = curve_to_q(tmp_beta);
    [q_tmp2,~,~,~] = Find_Rotation_and_Seed_unique(q_circ,q_tmp,1);
    d2(i) = InnerProd_Q(q_tmp2-q_circ,q_tmp2-q_circ);
    i
end

[~,idxc] = sort(d2,'ascend');

% Plot in order sequentially
for i=1:size(init_TOPI,3)
    clf
    imshow(I)
    figure(gcf)
    hold on
    plot(init_TOPI(2,:,idxc(i)),init_TOPI(1,:,idxc(i)),'r','LineWidth',5)
    pause
end

% Plot top 30 simultaneously
imshow(I)
figure(gcf)
hold on
for i=1:30
    plot(init_TOPI(2,:,idxc(i)),init_TOPI(1,:,idxc(i)),'r','LineWidth',5)
end

%% Run BAC using TOP initialization
n_curves = 30;
init_TOPIc = init_TOPI(:,:,idxc(1:n_curves));

% First estimates separate interior/exterior density for each contour,
% whereas second estimates a "pooled" interior/exterior density
[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegNT(I,init_TOPIc,0.1,1,0,500,[],[],0,0);
[seg{2},totE{2},interior{2},tmr{2}] = TOPBACSegNT(I,init_TOPIc,0.1,1,0,500,[],[],0,1);