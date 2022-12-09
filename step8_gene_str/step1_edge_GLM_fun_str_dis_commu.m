%% fitlm between str connectivity, dis, commu and gene co-expression
%% the dis and communicability are calculated by averaging individual auc
%% para
path = pwd;

roi_num = 48;
R_index = 8:2:48;

rand_time = 10000;
p_thr = 0.05;
p_fdr = 0.05;

results = struct;


%% data
% dis data
load('dis_whole.mat','dis_matrix_bi')

% node name
load('node_name.mat')

% gene data
[gene_data,gene_name,~] = xlsread([pwd '\gene_JHU.csv']);
gene_data(:,1) = [];
nan_region = find(isnan(gene_data(:,1)));
gene_name(1) = [];
gene_data = corr(gene_data');

% str data
cd(path)
cd ..
str_path = [pwd '\step0_network_construction\AHBA_Smooth_DBM'];

load([path '\Sub_bin_dis_str_DBM.mat']);
str_data = zeros(roi_num,roi_num,6);
str_commu = zeros(roi_num,roi_num,6);
str_dis = zeros(roi_num,roi_num,6);
for isub = [1 6]
    load([str_path '\JS_KSDENSITY_256_Signal_sub_' num2str(isub,'%04d') '.mat'])
    str_data(:,:,isub) = SIM;

    str_commu(:,:,isub) = sub_topo_bin.ac(:,:,isub);
    str_dis(:,:,isub) = sub_topo_bin.adis(:,:,isub);
end
for isub = 2:5
    load([str_path '\JS_KSDENSITY_256_Signal_sub_' num2str(isub,'%04d') '.mat'])
    SIM(R_index,:) = nan;
    SIM(:,R_index) = nan;
    str_data(:,:,isub) = SIM;

    sub_topo_bin.ac(R_index,:) = nan;
    sub_topo_bin.ac(:,R_index) = nan;
    sub_topo_bin.adis(R_index,:) = nan;
    sub_topo_bin.adis(:,R_index) = nan;
    str_commu(:,:,isub) = sub_topo_bin.ac(:,:,isub);
    str_dis(:,:,isub) = sub_topo_bin.adis(:,:,isub);
end
str_data = nanmean(str_data,3);
str_commu = nanmean(str_commu,3);
str_dis = nanmean(str_dis,3);


%% cal
roi_index = setdiff(1:roi_num,nan_region);
for iroi = 1:length(roi_index)
    disp(['Now calculating the data in region (' num2str(iroi) '/' num2str(length(roi_index)) ')   |' datestr(clock)])
    
    Y = gene_data(:,roi_index(iroi));
    Y([roi_index(iroi);nan_region]) = [];
    X = [str_data(:,roi_index(iroi)),str_commu(:,roi_index(iroi)),str_dis(:,roi_index(iroi))];
    X([roi_index(iroi);nan_region],:) = [];
    
    dis_moran = dis_matrix_bi;
    dis_moran(:,[roi_index(iroi);nan_region]) = [];
    dis_moran([roi_index(iroi);nan_region],:) = [];
    W = 1./dis_moran;
    W(1:length(W)+1:end) = 1;
    MEM = compute_mem(W);
    Y_moran = moran_randomization(Y,MEM,rand_time);
    Y_moran = squeeze(Y_moran);
    
    mdl = fitlm(X,Y);
    
    results(iroi).roi = roi_index(iroi);
    results(iroi).region_name = node_name{roi_index(iroi)};
    results(iroi).X = X;
    results(iroi).Y = Y;
    results(iroi).AR = mdl.Rsquared.Adjusted;
    
    results(iroi).AR_rand = zeros(rand_time,1);
    results(iroi).AR_moran = zeros(rand_time,1);
    for irand = 1:rand_time
        rand_order = randperm(length(Y));
        while isequal(rand_order,1:length(Y)), rand_order = randperm(length(Y)); end
        Y_rand = Y(rand_order);
        mdl = fitlm(X,Y_rand);
        results(iroi).AR_rand(irand) = mdl.Rsquared.Adjusted;
        
        mdl = fitlm(X,Y_moran(:,irand));
        results(iroi).AR_moran(irand) = mdl.Rsquared.Adjusted;
    end
    results(iroi).AR_p = (sum(results(iroi).AR_rand >= results(iroi).AR)+1)/(rand_time+1);
    results(iroi).AR_moran_p = (sum(results(iroi).AR_moran >= results(iroi).AR)+1)/(rand_time+1);
end

results(1).p_FDR = gretna_FDR(cat(1,results(:).AR_p),p_fdr);
if ~isempty(results(1).p_FDR)
    results(1).sig_region = cat(1,results(cat(1,results(:).AR_p) <= results(1).p_FDR).roi);
end

results(1).p_moran_FDR = gretna_FDR(cat(1,results(:).AR_moran_p),p_fdr);
if ~isempty(results(1).p_moran_FDR)
    results(1).sig_moran_region = cat(1,results(cat(1,results(:).AR_moran_p) <= results(1).p_moran_FDR).roi);
end


%% save
cd(path)
save('Edge_GLM_gene_str_dis_commu_DBM.mat','results')