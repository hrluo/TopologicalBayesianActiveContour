%%%%%%%%%%%%%%%
% Main Script %
%%%%%%%%%%%%%%%

%% Paths - start in D_home folder
% Replace D_home with your home path
D_home = '/Users/jdstrait/OneDrive - University of Georgia/Research/TOPBAC/COMPARE/REVISION/updated_code';
D_donut = 'Data/Donut';
D_lesion = 'Data/Lesion';
D_neuron = 'Data/Neuron';
D_bone = 'Data/Bone';
D_results = 'Results';

%% Simulated examples
% Donut example with Gaussian blur - Figures 2, 3
cd(D_donut)
load test_train.mat
load KarcherMeans.mat
cd(D_home)

init_TOPI{1} = 'donutblur_3_5_5.ppm';
init_TOPI{2} = 'donutblur_1_5_3.ppm';
init_TOPI{3} = 'donutblur_1_5_5.ppm';
init_TOPI{4} = 'donutblur_1_5_7.ppm';

[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegT(testIb,trainIb,train_beta,2,2,init_TOPI{1},0.15,0.3,0,[],1000,[],[],0.05);
[seg{2},totE{2},interior{2},tmr{2}] = TOPBACSegT(testIb,trainIb,train_beta,2,2,init_TOPI{2},0.15,0.3,0,[],1000,[],[],0.05);
[seg{3},totE{3},interior{3},tmr{3}] = TOPBACSegT(testIb,trainIb,train_beta,2,2,init_TOPI{3},0.15,0.3,0,[],1000,[],[],0.05);
[seg{4},totE{4},interior{4},tmr{4}] = TOPBACSegT(testIb,trainIb,train_beta,2,2,init_TOPI{4},0.15,0.3,0,[],1000,[],[],0.05);
[seg{5},totE{5},interior{5},tmr{5}] = TOPBACSegT(testIb,trainIb,train_beta,2,2,init_TOPI{1},0.15,0.3,0.05,[],1000,[],q_bar,0.05);
[seg{6},totE{6},interior{6},tmr{6}] = TOPBACSegT(testIb,trainIb,train_beta,2,2,init_TOPI{1},0.15,0.3,0.15,[],1000,[],q_bar,0.05);

% Appendix - bumpier initialization to show that aggressive image updates
% can form loops - Figure 13
[~,~,~,~,init_user] = TOPBACSegT(testI,trainI,train_beta,2,1,[],0.15,0.3,0,[],2);
for i=1:2
    init_TOPI{5}(:,:,i) = init_user{1};
end

[seg{7},totE{7},interior{7},tmr{7}] = TOPBACSegT(testI,trainI,train_beta,2,4,init_TOPI{5},0.15,0.3,0,[],1000);
[seg{8},totE{8},interior{8},tmr{8}] = TOPBACSegT(testI,trainI,train_beta,2,4,init_TOPI{5},0.4,0.3,0,[],1000);
[seg{9},totE{9},interior{9},tmr{9}] = TOPBACSegT(testI,trainI,train_beta,2,4,init_TOPI{5},0.8,0.3,0,[],1000);
[seg{10},totE{10},interior{10},tmr{10}] = TOPBACSegT(testI,trainI,train_beta,2,4,init_TOPI{5},0.4,3,0,[],1000);
[seg{11},totE{11},interior{11},tmr{11}] = TOPBACSegT(testI,trainI,train_beta,2,4,init_TOPI{5},0.8,3,0,[],1000);
[seg{12},totE{12},interior{12},tmr{12}] = TOPBACSegT(testI,trainI,train_beta,2,4,init_TOPI{5},0.8,10,0,[],1000);

% Evaluate all results using performance measures - Table 1
gtI = uint8(testI/255);
for i=1:12
    [HauD{i},JD{i},HamD{i},PM{i},ESD{i}] = SegDistTop(interior{i},gtI,seg{i},test_beta);
    i
end

cd(D_results)
save('Donut.mat','testIb','trainIb','testI','trainI','test_beta','train_beta','seg','totE','interior','tmr','HauD','JD','HamD','PM','ESD','init_TOPI')
clearvars -except D_home D_lesion D_results
cd(D_home)

%% Skin lesion examples
% Benign nevus example 1 (ISIC_0000424) - Figure 5
cd(D_lesion)
load BenignNevus1.mat
cd(D_home)
init_TOPI1 = 'ISIC_0000424_boundaries.ppm.bmp';
[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegT(testI,trainI,trainbeta,2,2,init_TOPI1,0.3,0.3,0,[],[],[],[],0.05);
cd(D_results)
save('BenignNevus1.mat','testI','seg','totE','interior','tmr')
clearvars -except D_home D_lesion D_results
cd(D_home)

% Benign nevus example 2 (ISIC_0000351) - Figure 6
cd(D_lesion)
load BenignNevus2.mat
cd(D_home)
init_TOPI2 = 'ISIC_0000351_boundaries.ppm.bmp';
[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegT(testI,trainI,trainbeta,2,2,init_TOPI2,0.3,0.3,0,[],[],[],[],0.05);
cd(D_results)
save('BenignNevus2.mat','testI','seg','totE','interior','tmr')
clearvars -except D_home D_lesion D_results
cd(D_home)

% Melanoma example 3 (ISIC_0000150) - Figure 14
cd(D_lesion)
load Melanoma.mat
cd(D_home)
init_TOPI3 = 'ISIC_0000150_boundaries.ppm.bmp';
[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegT(testI,trainI,trainbeta,3,2,init_TOPI3,[0.4,0.1,0.1],0.3,0,[],1000,[],[],0.05);
[seg{2},totE{2},interior{2},tmr{2}] = TOPBACSegT(testI,trainI,trainbeta,3,2,init_TOPI3,[0.4,0.1,0.1],0.3,0,[],1000,[],[],0.30);
cd(D_results)
save('Melanoma.mat','testI','seg','totE','interior','tmr')
clearvars -except D_home D_lesion D_results
cd(D_home)

% Benign nevus example 4 (ISIC_0000476) - Figure 15
cd(D_lesion)
load BenignNevus3.mat
cd(D_home)
init_TOPI4 = 'ISIC_0000476_boundaries.ppm.bmp';
[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegT(testI,trainI,trainbeta,2,2,init_TOPI4,0.3,0.3,0,[],1000,[],[],0.05);
[seg{2},totE{2},interior{2},tmr{2}] = TOPBACSegT(testI,trainI,trainbeta,2,2,init_TOPI4,0.3,0.3,0,[],1000,[],[],0.35);
cd(D_results)
save('BenignNevus3.mat','testI','seg','totE','interior','tmr')
clearvars -except D_home D_lesion D_results
cd(D_home)

%% Neuron example - Figure 7
% Note: code for pre-processing to obtain TOP initialization curves for the
% selected neuron image is found in the file NeuronScript.m
cd(D_neuron)
load NeuronPPD.mat
cd(D_home)

% First uses "pooled" density estimates, second uses individual density
% estimates with respect to each TOP initialized contour separately
n_curves = size(init_TOPIc,3);
lambda1 = [0.20*ones(1,25),0.10*ones(1,(n_curves-26+1))];
[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegNT(I,init_TOPIc,0.1,0.7,0,600,1e-8,0.30,1,0);
[seg{2},totE{2},interior{2},tmr{2}] = TOPBACSegNT(I,init_TOPIc,0.1,0.7,0,300,1e-8,0.05,1,0);
[seg{3},totE{3},interior{3},tmr{3}] = TOPBACSegNT(I,init_TOPIc,lambda1,0.7,0,1000,1e-8,0.30,1,0);
[seg{4},totE{4},interior{4},tmr{4}] = TOPBACSegNT(I,init_TOPIc,lambda1,0.7,0,1000,1e-8,0.05,1,0);
save('NeuronSeg.mat','I','init_TOPIc','TOPI','seg','totE','interior','tmr')
cd(D_home)

%% MPEG-7 bone examples - Figures 9, 10, 11, 12; Tables 2, 3
cd(D_bone)
load test_train.mat
load KarMeans.mat
cd(D_home)

init_TOPI = 'TOPI.ppm';
init_TOPIb = 'TOPIb.ppm';
init_TOPIp = 'TOPIp.ppm';

% Comparison of no noise, contour perturbation, Gaussian blur
% TOP+BAC
[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegT(testI,trainI,train_beta,1,2,init_TOPI,0.3,0.3,0);
[seg{2},totE{2},interior{2},tmr{2}] = TOPBACSegT(testIp,trainIp,train_betap,1,2,init_TOPIp,0.3,0.3,0);
[seg{3},totE{3},interior{3},tmr{3}] = TOPBACSegT(testIb,trainIb,train_beta,1,2,init_TOPIb,0.3,0.3,0,[],[],[],[],0.03);

% BAC with user-specified initialization
[~,~,~,~,init_user] = TOPBACSegT(testI,trainI,train_beta,1,1,[],0.3,0.3,0,[],2);

[seg{4},totE{4},interior{4},tmr{4}] = TOPBACSegT(testI,trainI,train_beta,1,4,init_user{1},0.3,0.3,0);
[seg{5},totE{5},interior{5},tmr{5}] = TOPBACSegT(testIp,trainIp,train_betap,1,4,init_user{1},0.3,0.3,0);
[seg{6},totE{6},interior{6},tmr{6}] = TOPBACSegT(testIb,trainIb,train_beta,1,4,init_user{1},0.3,0.3,0,[],[],[],[],0.03);

% TOP
[~,~,~,~,TOP{7}] = TOPBACSegT(testI,trainI,train_beta,1,2,init_TOPI,0.3,0.3,0,[],2);
[~,~,~,~,TOP{8}] = TOPBACSegT(testIp,trainIp,train_betap,1,2,init_TOPIp,0.3,0.3,0,[],2);
[~,~,~,~,TOP{9}] = TOPBACSegT(testIb,trainIb,train_beta,1,2,init_TOPIb,0.3,0.3,0,[],2);

for i=7:9
    seg{i} = TOP{i}{1};
    interior{i} = uint8(roipoly(testI,seg{i}(2,:),seg{i}(1,:)));
end

% Evaluate results using performance measures
gt_Ib = uint8(testI/255);
gt_Ibp = uint8(testIp/255);
idx1 = [1,3,4,6,7,9];
idx2 = [2,5,8];

for i=1:6
    [HauD{idx1(i)},JD{idx1(i)},HamD{idx1(i)},PM{idx1(i)},ESD{idx1(i)}] = SegDistTop(interior{idx1(i)},gt_Ib,seg{idx1(i)},test_beta);
    i
end

for i=1:3
    [HauD{idx2(i)},JD{idx2(i)},HamD{idx2(i)},PM{idx2(i)},ESD{idx2(i)}] = SegDistTop(interior{idx2(i)},gt_Ibp,seg{idx2(i)},test_beta);
    i
end   

% Demonstrate use of prior
[seg{10},totE{10},interior{10},tmr{10}] = TOPBACSegT(testIsp,trainIsp,train_beta,1,4,init_user{1},0.3,0.3,0.05,[],[],[],q_bar);
[seg{11},totE{11},interior{11},tmr{11}] = TOPBACSegT(testIp,trainIp,train_betap,1,4,init_user{1},0.3,0.3,0.05,[],[],[],qp_bar);

cd(D_results)
save('Bone.mat')
clearvars -except D_home D_lesion D_results
cd(D_home)