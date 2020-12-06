%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Donut counterexample for k-means %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths - start in D_home folder
D_home = '/Users/jdstrait/OneDrive - University of Georgia/Research/TOPBAC/COMPARE/REVISION/updated_code';
D_data = 'Data/Donut';

% Load original donut image data
cd(D_data)
load('test_train.mat')
cd(D_home)

I{1} = testI;
for k=1:9
    I{k+1} = trainI{k};
end

%% POSSIBLE EXAMPLE 1 (?)
for k=1:10
    % Binary masks for exterior and interior of donut
    ext_mask = (I{k}==0);
    int_mask = (I{k}==255);
    
    % Change value for exterior and interior (bring closer to each other)
    IG{k} = uint8(0.3*255*int_mask+0.7*255*ext_mask);
    
    % Patch to change pixel value of
    center = [640,512];
    radius = 50;
    for i=1:size(IG{k},1)
        for j=1:size(IG{k},2)
            if sqrt((i-center(1))^2+(j-center(2))^2)<radius
                IG{k}(i,j) = 0.5*255;
            end
        end
    end
    
    % Interesting - k-means with IB1{k} gives vastly different result to
    % IBG1{k} due to blurring, since the patch is equidistant between
    % interior and exterior pixel values originally, but blurring might
    % impact that
    IGB{k} = imgaussfilt(IG{k},0.5);
    
    [kIG{k},~] = imsegkmeans(IG{k},2);
    [kIGB{k},~] = imsegkmeans(IGB{k},2);
end

% Show k-means with 2 clusters result - completely different topology under
% blurring compared to original image!
k = 1;

imshow(255*(kIG{k}-1))
figure(gcf)

imshow(255*(kIGB{k}-1))
figure(gcf)

%% POSSIBLE EXAMPLE 2
for k=1:10
    % Binary masks for exterior and interior of donut
    ext_mask = (I{k}==0);
    int_mask = (I{k}==255);
    
    % Change value for exterior and interior (bring closer to each other)
    IG2{k} = uint8(0.3*255*int_mask+0.7*255*ext_mask);
    
    % Four patches to change pixel value of
    radius = 50;
    center(1,:) = [0.2*768,0.2*1024];
    center(2,:) = [0.2*768,0.8*1024];
    center(3,:) = [0.8*768,0.2*1024];
    center(4,:) = [0.8*768,0.8*1024];
    pix = linspace(0.35,0.65,4);
    for m=1:4
        for i=1:size(IG2{k},1)
            for j=1:size(IG2{k},2)
                if sqrt((i-center(m,1))^2+(j-center(m,2))^2)<radius
                    IG2{k}(i,j) = 255*pix(m);
                end
            end
        end
    end
    
    % Interesting - k-means with IB1{k} gives vastly different result to
    % IBG1{k} due to blurring, since the patch is equidistant between
    % interior and exterior pixel values originally, but blurring might
    % impact that
    [kIG2{k},~] = imsegkmeans(IG2{k},2);
    [kIG22{k},~] = imsegkmeans(IG2{k},3);
end

% Show k-means with 2 clusters result - completely different topology under
% blurring compared to original image!
k = 1;

imshow(255*(kIG2{k}-1))
figure(gcf)