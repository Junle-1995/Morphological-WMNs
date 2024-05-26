%% compare edge via TFNBS
load('group_data_DBM.mat')

roi_num = 48;

data_group{1,1} = data.HC;
data_group{2,1} = data.MS;
data_group{3,1} = data.NMO;

M = 10000;

Mask_net = ones(roi_num);
Mask_net(1:length(Mask_net)+1:end) = 0;

Pthr = 0.05;

E = 0.5; % extension enhancement parameters
H = 2.25; % height enhancement parameters

Cova = data.cov;

[TFNBS_Result] = gretna_TFNBS_ftest_MS_NMO(data_group, M, Mask_net, Pthr, E, H, Cova);

results = TFNBS_Result;
save('Edge_compare_DBM.mat','results')