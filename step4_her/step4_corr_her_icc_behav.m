%% corr between her and behav corr
%% para
roi_num = 48;
Mask_net = logical(tril(ones(roi_num),-1));

path = pwd;
rand_time = 10000;

% her data
load('heritability_matrix_p_VBM.mat','her')
data_her = her(Mask_net);    

% dis data
load('dis_whole.mat','dis_matrix_bi')
dis_moran = dis_matrix_bi;
W = 1./dis_moran;
W(1:length(W)+1:end) = 1;
MEM = compute_mem(W);

% moran her
her(1:length(her)+1:end) = 1;
data_her_moran = zeros(length(data_her),rand_time);
for irand = 1:rand_time
    her_moran = moran_randomization(her,MEM,1);
    data_her_moran(:,irand) = her_moran(Mask_net);
end

% behav
cd ..
cd('step3_behav_corr')
load('Edge_BH_PLS_VBM.mat')
data_behav = [results.weights,results(4).weights];
load('Edge_BH_prediction_PCA_80_VBM.mat')
data_behav = [data_behav,mean(results.Consensus_weights,2),mean(results(4).Consensus_weights,2)];
data_behav = abs(data_behav);

%% cal
results = struct;

% behav
results.name = 'behav';
results.rho = corr(data_her,data_behav,'type','Spearman');
results.rho_moran = corr(data_her_moran,data_behav,'type','Spearman');
results.moran_p = (sum(abs(results.rho_moran) >= repmat(abs(results.rho),rand_time,1))+1)/(rand_time+1);

%% save
cd(path)
save('her_icc_behav_corr_VBM.mat','results')