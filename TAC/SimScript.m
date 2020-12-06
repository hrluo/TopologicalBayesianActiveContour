%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulated image script %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths - start in D_home folder
D_home = '/Users/jdstrait/OneDrive - University of Georgia/Research/TOPBAC/COMPARE/REVISION/updated_code';
D_donut = 'Data/Donut';

%% Donut with no noise (Figure 5)
cd(D_donut)
load test_train.mat    % loads data
cd(D_home)
init_TOPI = 'TOPI_donut.ppm';  % TOP initialization

[seg{1},totE{1},interior{1},tmr{1}] = TOPBACSegT(testI,trainI,train_beta(:,:,:,1),1,1,[],0.3,0.3,0);
[seg{2},totE{2},interior{2},tmr{2}] = TOPBACSegT(testI,trainI,train_beta(:,:,:,2),1,1,init1,0.3,0.3,0);
[seg{3},totE{3},interior{3},tmr{3}] = TOPBACSegT(testI,trainI,train_beta,2,2,init_TOPI,0.3,0.3,0);

% Performance measure evaluation
[HauD,JD,HamD,PM,ESD] = SegDistTop(interior{3},testI,seg{3},test_beta);