%% para
path = pwd;

load('node_name')
load('Edge_compare_DBM.mat')
load('group_data_DBM.mat')

fdr = 0.05;
p_correct = 1/1001;

post = struct;

%% tidy
PT_FDR = gretna_FDR(results.TFNBS_t_pval_sig_edge + p_correct,fdr);
if isempty(PT_FDR), PT_FDR = -1; end
for i = 1:results.N_sig_edge
    post(i).data = {squeeze(data.HC(results.Reg_sig_edge(i,1),results.Reg_sig_edge(i,2),:));...
        squeeze(data.MS(results.Reg_sig_edge(i,1),results.Reg_sig_edge(i,2),:));...
        squeeze(data.NMO(results.Reg_sig_edge(i,1),results.Reg_sig_edge(i,2),:))};
    
    post(i).i_index = results.Reg_sig_edge(i,1);
    post(i).j_index = results.Reg_sig_edge(i,2);

    post(i).inode = node_name{results.Reg_sig_edge(i,1)};
    post(i).jnode = node_name{results.Reg_sig_edge(i,2)};
    
    post(i).F = results.TFNBS_f_score_sig_edge(i);
    post(i).PF = results.TFNBS_f_pval_sig_edge(i)+p_correct;
    if post(i).PF > fdr, error('?'), end
    
    post(i).T = results.TFNBS_t_score_sig_edge(i,:)';
    post(i).PT = results.TFNBS_t_pval_sig_edge(i,:)'+p_correct;
    
    sig_index = find(post(i).PT <= PT_FDR);
    post(i).sig_index = sig_index;
    
    if isequal(sig_index,1)
        if post(i).T(sig_index) > 0, post(i).pattern = 'HC > MS';
        else post(i).pattern = 'HC < MS'; end
    elseif isequal(sig_index,2)
        if post(i).T(sig_index) > 0, post(i).pattern = 'HC > NMO';
        else post(i).pattern = 'HC < NMO'; end
    elseif isequal(sig_index,3)
        if post(i).T(sig_index) > 0, post(i).pattern = 'MS > NMO';
        else post(i).pattern = 'MS < NMO'; end
    elseif isequal(sig_index,[1;2])
        if sum(post(i).T(sig_index) > 0) == 2, post(i).pattern = 'HC > MS,NMO';
        elseif sum(post(i).T(sig_index) < 0) == 2, post(i).pattern = 'HC < MS,NMO';
        else post(i).pattern = 'error?!';
        end
    elseif isequal(sig_index,[1;3])
        if post(i).T(1) > 0 && post(i).T(3) < 0, post(i).pattern = 'HC,NMO > MS';
        elseif post(i).T(1) < 0 && post(i).T(3) > 0, post(i).pattern = 'HC,NMO < MS';
        else post(i).pattern = 'error?!';
        end
    elseif isequal(sig_index,[2;3])
        if sum(post(i).T(sig_index) > 0) == 2, post(i).pattern = 'HC,MS > NMO';
        elseif sum(post(i).T(sig_index) < 0) == 2, post(i).pattern = 'HC,MS < NMO';
        else post(i).pattern = 'error?!';
        end
    elseif isequal(sig_index,[1;2;3])
        post(i).pattern = 'Differential';
    end
end

%% save
cd(path)
save('Edge_compare_tidy_DBM.mat','post')