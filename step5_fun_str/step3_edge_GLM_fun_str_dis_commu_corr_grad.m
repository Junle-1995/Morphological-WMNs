%% corr the coupling R2 or P-value with functional gradient
roi_num = 48;

load('Edge_GLM_fun_str_dis_commu_DBM.mat')
AR = cat(1,results(:).AR);
AR_p = cat(1,results(:).AR_moran_p);

load('unrelated_gradient.mat','Gt_fun')
fun_grad = Gt_fun.gradients{1, 1}(:,1);
data = fun_grad;
data_name = {'Grad1'};


rand_time = 10000;

load('dis_whole.mat', 'dis_matrix_bi')
dis = dis_matrix_bi;
W = 1./dis;
W(1:length(W)+1:end) = 1;
MEM = compute_mem(W);
fun_moran = moran_randomization(data,MEM,rand_time);

results = struct;

[r,p] = corr([AR,AR_p],data,'type','Spearman');
r_moran = zeros(size(r,1),size(r,2),rand_time);
for irand = 1:rand_time
    r_moran(:,:,irand) = corr([AR,AR_p],fun_moran(:,:,irand),'type','Spearman');
end

for i = 1:length(data_name)
    results(i).AR = 'AR';
    results(i).fun = data_name{i};
    results(i).AR_data = AR;
    results(i).fun_data = data(:,i);
    results(i).r = r(1,i);
    results(i).p = p(1,i);
    results(i).r_moran = squeeze(r_moran(1,i,:));
    results(i).p_moran = (sum(abs(results(i).r_moran) > abs(results(i).r))+1)/(rand_time+1);
    
    results(i+length(data_name)).AR = 'AR_p';
    results(i+length(data_name)).fun = data_name{i};
    results(i+length(data_name)).AR_data = AR_p;
    results(i+length(data_name)).fun_data = data(:,i);
    results(i+length(data_name)).r = r(2,i);
    results(i+length(data_name)).p = p(2,i);
    results(i+length(data_name)).r_moran = squeeze(r_moran(2,i,:));
    results(i+length(data_name)).p_moran = (sum(abs(results(i+length(data_name)).r_moran) > abs(results(i+length(data_name)).r))+1)/(rand_time+1);
end

save('Edge_GLM_fun_str_dis_commu_corr_grad_DBM.mat','results')