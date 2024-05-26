%% fitlm between str connectivity, dis, commu and chemoarchitectonic connectivity
%% the dis and communicability are calculated by averaging individual auc
%% para
roi_num = 48;

load('Sub_matrix_str_DBM.mat')
str_data = mean(sub_matrix_smo,3);
load('neuro_trans_net.mat')
fun_data = TM_net.JHU;
load('Sub_bin_dis_str_DBM.mat')
str_commu = mean(sub_topo_bin.ac,3);
str_dis = mean(sub_topo_bin.adis,3);
load('dis_whole.mat','dis_matrix_bi')

path = pwd;
cd ..
load('node_name.mat')

rand_time = 10000;
p_fdr = 0.05;

results = struct;


%% cal
for iroi = 1:roi_num
    disp(['Now calculating the data in region (' num2str(iroi) '/' num2str(roi_num) ')   |' datestr(clock)])
    
    Y = fun_data(:,iroi);
    Y(iroi) = [];
    X = [str_data(:,iroi),str_commu(:,iroi),str_dis(:,iroi)];
    X(iroi,:) = [];
    
    dis_moran = dis_matrix_bi;
    dis_moran(:,iroi) = [];
    dis_moran(iroi,:) = [];
    W = 1./dis_moran;
    W(1:length(W)+1:end) = 1;
    MEM = compute_mem(W);
    Y_moran = moran_randomization(Y,MEM,rand_time);
    Y_moran = squeeze(Y_moran);
    
    mdl = fitlm(X,Y);
    
    results(iroi).roi = iroi;
    results(iroi).region_name = node_name{iroi};
    results(iroi).X = X;
    results(iroi).Y = Y;
    results(iroi).AR = mdl.Rsquared.Adjusted;
    
    results(iroi).AR_rand = zeros(rand_time,1);
    results(iroi).AR_moran = zeros(rand_time,1);
    for irand = 1:rand_time
        rand_order = randperm(roi_num-1);
        while isequal(rand_order,1:roi_num-1), rand_order = randperm(roi_num-1); end
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
    results(1).sig_region = find(cat(1,results(:).AR_p) <= results(1).p_FDR);
end

results(1).p_moran_FDR = gretna_FDR(cat(1,results(:).AR_moran_p),p_fdr);
if ~isempty(results(1).p_moran_FDR)
    results(1).sig_moran_region = find(cat(1,results(:).AR_moran_p) <= results(1).p_moran_FDR);
end


%% save
cd(path)
save('Edge_GLM_fun_str_dis_commu_DBM.mat','results')