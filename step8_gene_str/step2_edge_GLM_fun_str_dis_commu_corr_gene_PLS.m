%% corr the coupling R2 or P-value with gene expression via PLS
%% para
path = pwd;

pls_dim = 1;
roi_num = 48;
rand_time = 10000;
p_thr = 0.05;

gen_results = struct;


%% data
% AR data
load('Edge_GLM_gene_str_dis_commu_DBM.mat')
roi_index = cat(1,results(:).roi);
AR = cat(1,results(:).AR);
AR_p = -log10(cat(1,results(:).AR_moran_p));

% gene data
[gene_data,~,~] = xlsread([pwd '\gene_JHU.csv']);
gene_data(:,1) = [];
nan_region = find(isnan(gene_data(:,1)));
if ~isequal(nan_region,setdiff(1:roi_num,roi_index)')
    error('?')
end
gene_data(nan_region,:) = [];

% dis data for moran
load('dis_whole.mat','dis_matrix_bi')
dis = dis_matrix_bi;
dis(nan_region,:) = [];
dis(:,nan_region) = [];
W = 1./dis;
W(1:length(W)+1:end) = 1;
MEM = compute_mem(W);

%% PLS
gen_results(1).type = 'AR';
gen_results(2).type = 'AR_p';

% real
[~,~,XS,YS,~,PCTVAR,~,stats] = plsregress(zscore(gene_data),AR,pls_dim);
gen_results(1).data = AR;
gen_results(1).YS = YS;
gen_results(1).PCTVAR = PCTVAR;
gen_results(1).X_Score = XS;
gen_results(1).gen_weight = stats.W;

[~,~,XS,YS,~,PCTVAR,~,stats] = plsregress(zscore(gene_data),AR_p,pls_dim);
gen_results(2).data = AR_p;
gen_results(2).YS = YS;
gen_results(2).PCTVAR = PCTVAR;
gen_results(2).X_Score = XS;
gen_results(2).gen_weight = stats.W;

% moran
gen_results(1).PCTVAR_rand = zeros(rand_time,pls_dim);
gen_results(2).PCTVAR_rand = zeros(rand_time,pls_dim);
for irand = 1:rand_time
    disp(['permutation ' num2str(irand)...
        ' of ' num2str(rand_time) ' in AR   |' datestr(clock)]);

    gene_data_moran = moran_randomization(gene_data,MEM,1);

    [~,~,~,~,~,PCTVAR,~,stats] = plsregress(zscore(gene_data_moran),AR,pls_dim);
    gen_results(1).PCTVAR_rand(irand,:) = PCTVAR(2,:);
    gen_results(1).gen_weight_moran(:,irand) = stats.W;

    [~,~,~,~,~,PCTVAR,~,stats] = plsregress(zscore(gene_data_moran),AR_p,pls_dim);
    gen_results(2).PCTVAR_rand(irand,:) = PCTVAR(2,:);
    gen_results(2).gen_weight_moran(:,irand) = stats.W;
end
gen_results(1).p_rand = ...
    (sum(gen_results(1).PCTVAR_rand >= repmat(gen_results(1).PCTVAR(2,:),rand_time,1))+1)/(rand_time + 1);
gen_results(2).p_rand = ...
    (sum(gen_results(2).PCTVAR_rand >= repmat(gen_results(2).PCTVAR(2,:),rand_time,1))+1)/(rand_time + 1);


%% save
cd(path)
gen_results1 = gen_results;
gen_results = gen_results1(1);
save('Edge_GLM_gene_str_dis_commu_gene_PLS_DBM.mat','gen_results')
gen_results = gen_results1(2);
save('Edge_GLM_gene_str_dis_commu_gene_PLS_DBM_P.mat','gen_results')