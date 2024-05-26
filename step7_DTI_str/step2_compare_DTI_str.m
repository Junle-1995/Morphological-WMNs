%% diff between str connectivity with DTI fiber number (> 0 and =< 0)
%% para
path = pwd;

load('unrelated_sub_matrix_DTI.mat')
DTI_data = DTI_matrix;

Type = {'DBM','VBM'};

cd ..
load('node_name.mat')
roi_num = 48;

DTI_thr = 0;
p_fdr = 0.05;

results = struct;


%% cal
for itype = 1:length(Type)
    load(['unrelated_sub_matrix_str_' Type{itype} '.mat'])
    str_data = sub_matrix;

    results(itype).Type = Type{itype};
    for isub = 1:size(str_data,3)
        DTI_i = DTI_data(:,:,isub);
        str_i = str_data(:,:,isub);

        DTI_i(1:length(DTI_i)+1:end) = nan;
        results(itype).data(isub,:) = [mean(str_i(DTI_i <= DTI_thr)),mean(str_i(DTI_i > DTI_thr))];
    end
    [h,p,ci,stats] = ttest(results(itype).data(:,1), results(itype).data(:,2));
    results(itype).T = stats.tstat;
    results(itype).P = p;
end


%% save
cd(path)
save('results_compare_str_DTI.mat','results')