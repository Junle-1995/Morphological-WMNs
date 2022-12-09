%% spatial corr between smo and no_smo
load('Sub_matrix_str_DBM.mat')
sub_matrix_no_smo_DBM = mean(sub_matrix_no_smo,3);
sub_matrix_smo_DBM = mean(sub_matrix_smo,3);

load('Sub_matrix_str_VBM.mat')
sub_matrix_no_smo_VBM = mean(sub_matrix_no_smo,3);
sub_matrix_smo_VBM = mean(sub_matrix_smo,3);

load('dis_whole.mat', 'dis_matrix_bi')

net_mask = logical(tril(ones(length(sub_matrix_no_smo_DBM)),-1));

rand_time = 10000;

dis = dis_matrix_bi;
W = 1./dis;
W(1:length(W)+1:end) = 1;
MEM = compute_mem(W);


%% corr
% DBM_smo_no_smo
sub_matrix_moran = moran_randomization(sub_matrix_no_smo_DBM,MEM,rand_time);

results(1).name = 'DBM_smo_no_smo';
results(1).no_smo = sub_matrix_no_smo_DBM;
results(1).smo = sub_matrix_smo_DBM;
results(1).mask = net_mask;
[results(1).r,results(1).p] = corr(sub_matrix_no_smo_DBM(net_mask),sub_matrix_smo_DBM(net_mask),'type','Spearman');

for irand = 1:rand_time
    data_moran = sub_matrix_moran(:,:,irand);
    results(1).r_rand(irand,1) = corr(data_moran(net_mask),sub_matrix_smo_DBM(net_mask),'type','Spearman');
end
results(1).p_rand = (sum(abs(results(1).r_rand) >= abs(results(1).r))+1)/(rand_time+1);

% VBM_smo_no_smo
sub_matrix_moran = moran_randomization(sub_matrix_no_smo_VBM,MEM,rand_time);

results(2).name = 'VBM_smo_no_smo';
results(2).no_smo = sub_matrix_no_smo_VBM;
results(2).smo = sub_matrix_smo_VBM;
results(2).mask = net_mask;
[results(2).r,results(2).p] = corr(sub_matrix_no_smo_VBM(net_mask),sub_matrix_smo_VBM(net_mask),'type','Spearman');

for irand = 1:rand_time
    data_moran = sub_matrix_moran(:,:,irand);
    results(2).r_rand(irand,1) = corr(data_moran(net_mask),sub_matrix_smo_VBM(net_mask),'type','Spearman');
end
results(2).p_rand = (sum(abs(results(2).r_rand) >= abs(results(2).r))+1)/(rand_time+1);

% smo_DBM_VBM
sub_matrix_moran = moran_randomization(sub_matrix_smo_DBM,MEM,rand_time);

results(3).name = 'smo_DBM_VBM';
results(3).no_smo = sub_matrix_smo_DBM;
results(3).smo = sub_matrix_smo_VBM;
results(3).mask = net_mask;
[results(3).r,results(3).p] = corr(sub_matrix_smo_DBM(net_mask),sub_matrix_smo_VBM(net_mask),'type','Spearman');

for irand = 1:rand_time
    data_moran = sub_matrix_moran(:,:,irand);
    results(3).r_rand(irand,1) = corr(data_moran(net_mask),sub_matrix_smo_VBM(net_mask),'type','Spearman');
end
results(3).p_rand = (sum(abs(results(3).r_rand) >= abs(results(3).r))+1)/(rand_time+1);

% no_smo_DBM_VBM
sub_matrix_moran = moran_randomization(sub_matrix_no_smo_DBM,MEM,rand_time);

results(4).name = 'no_smo_DBM_VBM';
results(4).no_smo = sub_matrix_no_smo_DBM;
results(4).smo = sub_matrix_no_smo_VBM;
results(4).mask = net_mask;
[results(4).r,results(4).p] = corr(sub_matrix_no_smo_DBM(net_mask),sub_matrix_no_smo_VBM(net_mask),'type','Spearman');

for irand = 1:rand_time
    data_moran = sub_matrix_moran(:,:,irand);
    results(4).r_rand(irand,1) = corr(data_moran(net_mask),sub_matrix_no_smo_VBM(net_mask),'type','Spearman');
end
results(4).p_rand = (sum(abs(results(4).r_rand) >= abs(results(4).r))+1)/(rand_time+1);

save('Edge_corr.mat','results')