%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate donut images %
%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths - start in D_home folder
D_home = '/Users/jdstrait/OneDrive - University of Georgia/Research/TOPBAC/COMPARE/REVISION/updated_code';
D_donut = 'Data/Donut';
D_source = 'Data/Donut/Original_Images';

%% Load in source images
set(0,'DefaultFigureVisible','off')
cd(D_source)
S = dir(fullfile([D_home,'/',D_source],'donut-*.ppm')); % pattern to match filenames
N = 200;  % number of points to re-sample curves to

%% Set up data
for k=1:numel(S)
    % Start with original MPEG-7 test image, and extract the boundary
    F = fullfile([D_home,'/',D_source],S(k).name);
    [filepath,fname,ext] = fileparts(S(k).name);
    fname = string(fname);
    origI = imread(F);
    origI = im2double(origI(:,:,1));
    tmp = bwboundaries(origI);
    
    orig_bdy1 = tmp{1}';
    orig_bdy2 = tmp{2}';
    cd(D_home)
    orig_bdyRS(:,:,1) = ReSampleCurve(orig_bdy1,N);
    orig_bdyRS(:,:,2) = ReSampleCurve(orig_bdy2,N);
    cd(D_donut)
    
    % Use original extracted boundary to create modified image (this is
    % done so that the image has more black pixels surrounding the target
    % object, as well as to be consistent with how the images with Gaussian
    % contour perturbation are generated), and extract true boundary
    fig = figure;
    hold on
    for i=1:2
        plot(orig_bdyRS(1,:,i),orig_bdyRS(2,:,i),'b','LineWidth',2);
    end
    axis off equal;
    f = getframe(fig);
    img = f.cdata;
    img = imresize(img,[768 NaN]);
    img = rgb2gray(img);
    img = imbinarize(img);
    img_out = imfill(~img,'holes');
    img_in = imfill(~img,[size(img,1)/2 size(img,2)/2]);
    img = abs(img_in-img_out);
    I{k} = uint8(255*img);
    tmp = bwboundaries(I{k});
    
    cd(D_home)
    for i=1:2
        true_bdy{k}(:,:,i) = ReSampleCurve(tmp{i}',N);
    end
    
    % Salt and pepper noise images
    Isp{k} = imnoise(I{k},'salt & pepper',0.3);
    
    % Gaussian blur images
    Ib{k} = imgaussfilt(I{k},15);
    
    % Gaussian contour perturbation images
    for i=1:2
        pert_bdyRS(:,:,i) = orig_bdyRS(:,:,i) + normrnd(0,3,[2,N]);
        pert_bdyRS(:,end,i) = pert_bdyRS(:,1,i);
    end
    fig = figure;
    hold on
    for i=1:2
        plot(pert_bdyRS(1,:,i),pert_bdyRS(2,:,i),'LineWidth',2);
    end
    axis off equal;
    f = getframe(fig);
    img = f.cdata;
    img = imresize(img,[768 NaN]);
    img = rgb2gray(img);
    img = imbinarize(img);
    img_out = imfill(~img,'holes');
    img_in = imfill(~img,[size(img,1)/2 size(img,2)/2]);
    img = abs(img_in-img_out);
    Ip{k} = uint8(255*img);
    
    tmp = bwboundaries(Ip{k});
    for i=1:length(tmp)
        rng = range(tmp{i});
        ArI(i) = rng(1)*rng(2);
    end
    [~,idx] = sort(ArI,'descend');
    clear ArI
    
    for i=1:2
        true_bdyp{k}(:,:,i) = ReSampleCurve(tmp{idx(i)}',N);
    end
    
    k
end

% Separate into test image versus training image set (and equivalent
% subsets for curves, into format easy for reading into TOPBACSeg program),
% and save relevant files (.mat and .ppm for TOP initialization)
idx_test = 1; 
idx_train = 1:numel(S); idx_train(idx_test) = [];

testI = I{idx_test}; test_beta = true_bdy{idx_test};
testIsp = Isp{idx_test};
testIp = Ip{idx_test}; test_betap = true_bdyp{idx_test};
testIb = Ib{idx_test};

trainI = I(idx_train);
trainIsp = Isp(idx_train);
trainIp = Ip(idx_train);
trainIb = Ib(idx_train);
for i=1:2
    for j=1:length(idx_train)
        train_beta(:,:,j,i) = true_bdy{idx_train(j)}(:,:,i);
        train_betap(:,:,j,i) = true_bdyp{idx_train(j)}(:,:,i);
    end
end

cd(D_donut)
save('test_train.mat','testI','testIsp','testIp','testIb','trainI','trainIsp','trainIp','trainIb','test_beta','test_betap','train_beta','train_betap')
imwrite(testI,'orig_test.ppm')
imwrite(testIsp,'SP_test.ppm')
imwrite(testIp,'GP_test.ppm')
imwrite(testIb,'GB_test.ppm')
cd(D_home)

%% k-means initialization/segmentation
[kI,~] = imsegkmeans(testI,2);
[kIsp,~] = imsegkmeans(testIsp,2);
[kIp,~] = imsegkmeans(testIp,2);
[kIb,~] = imsegkmeans(testIb,2);

% Obtain binary mask of clusters (for evaluation)
kI_mask = kI-1;
kIsp_mask = kIsp-1;
kIp_mask = kIp-1;
kIb_mask = kIb-1;

% Obtain a boundary estimate (for use with initialization), using same
% criterion as in TOPBACSeg.m program (crude estimate of area for each
% detected boundary, choose the largest)
tmp = bwboundaries(kI_mask);
new_boundary = tmp;
for i=1:length(tmp)
    rng = range(new_boundary{i});
    ArI(i) = rng(1)*rng(2);
end
[~,idx] = sort(ArI,'descend');
for i=1:2
    init_kI(:,:,i) = ReSampleCurve(new_boundary{idx(i)}',N);
end

tmp = bwboundaries(kIsp_mask);
new_boundary = tmp;
for i=1:length(tmp)
    rng = range(new_boundary{i});
    ArIsp(i) = rng(1)*rng(2);
end
[~,idx] = sort(ArIsp,'descend');
for i=1:2
    init_kIsp(:,:,i) = ReSampleCurve(new_boundary{idx(i)}',N);
end

tmp = bwboundaries(kIp_mask);
new_boundary = tmp;
for i=1:length(tmp)
    rng = range(new_boundary{i});
    ArIp(i) = rng(1)*rng(2);
end
[~,idx] = sort(ArIp,'descend');
for i=1:2
    init_kIp(:,:,i) = ReSampleCurve(new_boundary{idx(i)}',N);
end

tmp = bwboundaries(kIb_mask);
new_boundary = tmp;
for i=1:length(tmp)
    rng = range(new_boundary{i});
    ArIb(i) = rng(1)*rng(2);
end
[~,idx] = sort(ArIb,'descend');
for i=1:2
    init_kIb(:,:,i) = ReSampleCurve(new_boundary{idx(i)}',N);
end

cd(D_data)
save('K_means.mat','kI_mask','kIsp_mask','kIp_mask','kIb_mask','init_kI','init_kIsp','init_kIp','init_kIb')
cd(D_home)


%% below is un-verified
%% TOP initialization/segmentation
fI = fullfile(D_img_mod,'orig_test.ppm');
command1 = './mean_shift_v6/mean_shift ./'+string(fI)+' 1 5 5'
system(command1)
pause(5)
command2 = 'mv ./boundaries.ppm ./TOPI.ppm'
system(command2)
pause(5)
TOPI_mask = imread('./TOPI.ppm');
TOPI_mask = TOPI_mask(:,:,1);

fIsp = fullfile(D_img_mod,'SP_test.ppm');
command1 = './mean_shift_v6/mean_shift ./'+string(fIsp)+' 1 5 5'
system(command1)
pause(5)
command2 = 'mv ./boundaries.ppm ./TOPIsp.ppm'
system(command2)
pause(5)
TOPIsp_mask = imread('./TOPIsp.ppm');
TOPIsp_mask = TOPIsp_mask(:,:,1);

fIp = fullfile(D_img_mod,'GP_test.ppm');
command1 = './mean_shift_v6/mean_shift ./'+string(fIp)+' 1 5 5'
system(command1)
pause(5)
command2 = 'mv ./boundaries.ppm ./TOPIp.ppm'
system(command2)
pause(5)
TOPIp_mask = imread('./TOPIp.ppm');
TOPIp_mask = TOPIp_mask(:,:,1);

fIb = fullfile(D_img_mod,'GB_test.ppm');
command1 = './mean_shift_v6/mean_shift ./'+string(fIb)+' 1 5 5'
system(command1)
pause(5)
command2 = 'mv ./boundaries.ppm ./TOPIb.ppm'
system(command2)
pause(5)
TOPIb_mask = imread('./TOPIb.ppm');
TOPIb_mask = TOPIb_mask(:,:,1);

%% Compute elastic shape means of true training image boundaries
cd(D_BAC)
for i=1:2
    [~,q_bar(:,:,i)] = FindElasticMean(train_beta(:,:,:,i),0);
    [~,qp_bar(:,:,i)] = FindElasticMean(train_betap(:,:,:,i),0);
end
cd(D_img_mod)
save('KarMeans.mat','q_bar','qp_bar')
cd(D_home)