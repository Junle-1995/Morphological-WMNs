%% corr the coupling R2 or P-value with transmitter maps
%% para
roi_num = 48;

load('Edge_GLM_fun_str_dis_commu_DBM.mat')
AR = cat(1,results(:).AR);
AR_p = cat(1,results(:).AR_moran_p);

load('neuro_trans_net.mat','TM_data','trans_name')
TM_data = TM_data.JHU;

p_fdr = 0.05;
rand_time = 10000;

results_trans = struct;

%% moran
load('dis_whole.mat', 'dis_matrix_bi')
dis = dis_matrix_bi;
W = 1./dis;
W(1:length(W)+1:end) = 1;

trans_moran = nan(size(TM_data,1),size(TM_data,2),rand_time);
for itrans = 1:size(TM_data,2)
    TM_data_i = TM_data(:,itrans);
    nan_region = isnan(TM_data_i);
    TM_data_i(nan_region) = [];
    
    W1 = W;
    W1(nan_region,:) = [];
    W1(:,nan_region) = [];
    MEM = compute_mem(W1);
    trans_moran(~nan_region,itrans,:) = moran_randomization(TM_data_i,MEM,rand_time);
end


%% cal
[r_trans,p_trans] = corr([AR,AR_p],TM_data,'type','Spearman','rows','pairwise');
r_trans_moran = zeros(size(r_trans,1),size(r_trans,2),rand_time);
for irand = 1:rand_time
    r_trans_moran(:,:,irand) = corr([AR,AR_p],trans_moran(:,:,irand),'type','Spearman','rows','pairwise');
end

results_trans(1).AR = 'AR';
results_trans(1).AR_data = AR;
results_trans(1).fun_data = TM_data;
results_trans(1).r = r_trans(1,:)';
results_trans(1).p = p_trans(1,:)';
results_trans(1).r_moran = squeeze(r_trans_moran(1,:,:));
results_trans(1).p_moran = (sum(abs(results_trans(1).r_moran) > repmat(abs(results_trans(1).r),1,rand_time),2)+1)/(rand_time+1);
results_trans(1).p_moran_FDR = gretna_FDR(results_trans(1).p_moran,p_fdr);
if ~isempty(results_trans(1).p_moran_FDR)
    results_trans(1).sig_trans = trans_name(results_trans(1).p_moran <= results_trans(1).p_moran_FDR);
    results_trans(1).sig_trans(:,2) = num2cell(results_trans(1).r(results_trans(1).p_moran <= results_trans(1).p_moran_FDR));
end

results_trans(2).AR = 'AR_p';
results_trans(2).AR_data = AR_p;
results_trans(2).fun_data = TM_data;
results_trans(2).r = r_trans(2,:)';
results_trans(2).p = p_trans(2,:)';
results_trans(2).r_moran = squeeze(r_trans_moran(2,:,:));
results_trans(2).p_moran = (sum(abs(results_trans(2).r_moran) > repmat(abs(results_trans(2).r),1,rand_time),2)+1)/(rand_time+1);
results_trans(2).p_moran_FDR = gretna_FDR(results_trans(2).p_moran,p_fdr);
if ~isempty(results_trans(2).p_moran_FDR)
    results_trans(2).sig_trans = trans_name(results_trans(2).p_moran <= results_trans(2).p_moran_FDR);
    results_trans(2).sig_trans(:,2) = num2cell(results_trans(2).r(results_trans(2).p_moran <= results_trans(2).p_moran_FDR));
end

save('Edge_GLM_fun_str_dis_commu_corr_trans_DBM.mat','results_trans','results_degree')