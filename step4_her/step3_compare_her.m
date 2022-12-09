%% compare the heritability between DBM and VBM
path = pwd;

roi_num = 48;
her_mask = logical(tril(ones(roi_num),-1));

dbm_her_file = 'heritability_matrix_p_DBM.mat';
vbm_her_file = 'heritability_matrix_p_VBM.mat';

load(dbm_her_file)
dbm_her = her(her_mask);
dbm_p = her_p(her_mask);

load(vbm_her_file)
vbm_her = her(her_mask);
vbm_p = her_p(her_mask);

cd(path)
[~,p,~,stats] = ttest(dbm_her,vbm_her);