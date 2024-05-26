%% corr the duration and EDSS with altered edge, controling the cov in the edge
%% para
load('Edge_compare_tidy_DBM.mat')
load('group_data_DBM.mat')

MS_cov = data.cov{2};
MS_cov(:,1) = tiedrank(MS_cov(:,1));
MS_cov(:,2) = tiedrank(MS_cov(:,2));
NMO_cov = data.cov{3};
NMO_cov(:,1) = tiedrank(NMO_cov(:,1));
NMO_cov(:,2) = tiedrank(NMO_cov(:,2));

MS_duration = tiedrank(data.cli{2}(:,1));
MS_EDSS = tiedrank(data.cli{2}(:,2));
NMO_duration = tiedrank(data.cli{3}(:,1));
NMO_EDSS = tiedrank(data.cli{3}(:,2));

df_MS = size(MS_cov,1)-2-rank(MS_cov);
df_NMO = size(NMO_cov,1)-2-rank(NMO_cov);

fdr = 0.05;


%% bin
MS_corr = struct;
NMO_corr = struct;
for i = 1:length(post)    
    % MS
    MS_corr(i).i_index = post(i).i_index;
    MS_corr(i).j_index = post(i).j_index;
    MS_corr(i).inode = post(i).inode;
    MS_corr(i).jnode = post(i).jnode;
    MS_corr(i).data = post(i).data{2};
    bi = regress(MS_corr(i).data,[ones(size(MS_corr(i).data,1),1) MS_cov]);
    MS_corr(i).data_res = MS_corr(i).data - MS_cov*bi(2:end);
    
    if sum(post(i).sig_index == 1) == 1
        % duration
        [MS_corr(i).rho_duration,~] = corr(MS_corr(i).data_res,MS_duration);
        t = abs(MS_corr(i).rho_duration*sqrt(df_MS/(1-MS_corr(i).rho_duration^2)));
        MS_corr(i).p_duration = 2*(1-tcdf(t,df_MS)); % two-tail
        
        % EDSS
        [MS_corr(i).rho_EDSS,~] = corr(MS_corr(i).data_res,MS_EDSS);
        t = abs(MS_corr(i).rho_EDSS*sqrt(df_MS/(1-MS_corr(i).rho_EDSS^2)));
        MS_corr(i).p_EDSS = 2*(1-tcdf(t,df_MS)); % two-tail
    end
    
    
    % NMO
    NMO_corr(i).inode = post(i).inode;
    NMO_corr(i).jnode = post(i).jnode;
    NMO_corr(i).data = post(i).data{3};
    bi = regress(NMO_corr(i).data,[ones(size(NMO_corr(i).data,1),1) NMO_cov]);
    NMO_corr(i).data_res = NMO_corr(i).data - NMO_cov*bi(2:end);
    
    if sum(post(i).sig_index == 2) == 1
        % duration
        [NMO_corr(i).rho_duration,~] = corr(NMO_corr(i).data_res,NMO_duration);
        t = abs(NMO_corr(i).rho_duration*sqrt(df_NMO/(1-NMO_corr(i).rho_duration^2)));
        NMO_corr(i).p_duration = 2*(1-tcdf(t,df_NMO)); % two-tail
        
        % EDSS
        [NMO_corr(i).rho_EDSS,~] = corr(NMO_corr(i).data_res,NMO_EDSS);
        t = abs(NMO_corr(i).rho_EDSS*sqrt(df_NMO/(1-NMO_corr(i).rho_EDSS^2)));
        NMO_corr(i).p_EDSS = 2*(1-tcdf(t,df_NMO)); % two-tail
    end
end
MS_p_all = [cat(1,MS_corr(:).p_duration);cat(1,MS_corr(:).p_EDSS)];
NMO_p_all = [cat(1,NMO_corr(:).p_duration);cat(1,NMO_corr(:).p_EDSS)];

MS_corr(1).p_FDR = gretna_FDR(MS_p_all,fdr);
NMO_corr(1).p_FDR = gretna_FDR(NMO_p_all,fdr);


%% save
save('Edge_corr_cli_DBM.mat','MS_corr','NMO_corr')