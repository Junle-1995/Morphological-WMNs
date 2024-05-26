%% generate the nii data of connectivity for two-way ranova
%% each nii contain data:3(icc type)
%% do two-way ranova in MRM toolbox
%% the file sort in MRM is: smo_dbm, smo_vbm, no_smo_dbm, no_smo_vbm
%% do the post hoc with paired t test
%% para
load('Sub_matrix_str_DBM.mat')
data.no_smo_dbm = mean(sub_matrix_no_smo,3);
data.smo_dbm = mean(sub_matrix_smo,3);

load('Sub_matrix_str_VBM.mat')
data.no_smo_vbm = mean(sub_matrix_no_smo,3);
data.smo_vbm = mean(sub_matrix_smo,3);

type_rank = {'smo_dbm','smo_vbm','no_smo_dbm','no_smo_vbm'};

roi_num = 48;
net_mask = logical(tril(ones(roi_num),-1));
[index_con_i,index_con_j] = find(net_mask);

contrast_name = {'Net_type','Smo_type','Smo_type_x_Net_type'};

p_thr = 0.05;

path = pwd;

results = struct;

post = struct;

%% generate
cd(path)
cd ..
cd('step0_network_construction\JHU')
atlas_path = pwd;

cd(path)
save_path = [path '\data_connectivity_ranova'];
mkdir(save_path)

[Vout, Y_mode, ~] = gretna_read_image([atlas_path '\JHU-WhiteMatter-labels-3mm.nii']);
Vout.dt(1) = 16;
for i = 1:length(type_rank)
    mkdir(save_path,type_rank{i})
    
    for icon = 1:sum(net_mask(:))
        Y_data = Y_mode(:,:,1:2);
        Y_data(:) = 0;
        Vdata = Vout;
        Vdata.dt(1) = 16;
        Vdata.dim(3) = 2;
        
        Y_data(1) = data.(type_rank{i})(index_con_i(icon),index_con_j(icon));
        
        gretna_write_image(Y_data, Vdata, [save_path '\' type_rank{i} '\' num2str(icon,'%04d') '.nii'])
    end
end
cd(path)

%% analysis the data in MRM
load('MRM_con.mat')
mkdir('data_connectivity_ranova_results')
MRM_estimate

%% collect the results
MRM_path = [path '\data_connectivity_ranova_results\Contrasts'];

for i = 1:length(contrast_name)
    cd([MRM_path '\' contrast_name{i}])

    % F
    [~, Ydata, ~] = gretna_read_image(['MRM_' contrast_name{i} '_F.nii']);
    results.(['F_' contrast_name{i}]) = Ydata(1);

    % P
    [~, Ydata, ~] = gretna_read_image(['MRM_' contrast_name{i} '_P_UC_.nii']);
    results.(['P_' contrast_name{i}]) = Ydata(1);

    % DF
    results.(['DF_' contrast_name{i}]) = load(['MRM_' contrast_name{i} '_DF.txt']);
end
cd(path)

%% post hoc
for i = 1:length(results)
    % Net type
    if results(i).P_Net_type <= p_thr
        data_dbm = [data.no_smo_dbm(net_mask);data.smo_dbm(net_mask)];
        
        data_vbm = [data.no_smo_vbm(net_mask);data.smo_vbm(net_mask)];
        
        [~,post(i).Net_p,~,stats] = ttest(data_dbm,data_vbm);
        post(i).Net_t = stats.tstat;
        post(i).Net_pattern = '';
        if post(i).Net_p < p_thr
            if post(i).Net_t > 0, post(i).Net_pattern = 'DBM > VBM';
            else post(i).Net_pattern = 'DBM < VBM'; end
        end
    end
    
    % Smo type
    if results(i).P_Smo_type <= p_thr
        data_smo = [data.smo_dbm(net_mask);data.smo_vbm(net_mask)];
        
        data_no_smo = [data.no_smo_dbm(net_mask);data.no_smo_vbm(net_mask)];
        
        [~,post(i).Smo_p,~,stats] = ttest(data_smo,data_no_smo);
        post(i).Smo_t = stats.tstat;
        post(i).Smo_pattern = '';
        if post(i).Smo_p < p_thr
            if post(i).Smo_t > 0, post(i).Smo_pattern = 'Smo > No_smo';
            else post(i).Smo_pattern = 'Smo < No_smo'; end
        end
    end
    
    
    % interaction
    if results(i).P_Smo_type_x_Net_type <= p_thr
        data_smo_dbm = data.smo_dbm(net_mask);
        data_smo_vbm = data.smo_vbm(net_mask);
        
        data_no_smo_dbm = data.no_smo_dbm(net_mask);
        data_no_smo_vbm = data.no_smo_vbm(net_mask);
        
        post(i).Inter_pattern = '';
        post(i).Inter_test = 'Smo:dbm_vbm; No_smo:dbm_vbm; DBM:smo_no_smo; VBM:smo_no_smo';
        
        % smo: dbm_vbm
        [~,post(i).Inter_p(1),~,stats] = ttest(data_smo_dbm,data_smo_vbm);
        post(i).Inter_t(1) = stats.tstat;
        
        % no_smo: dbm_vbm
        [~,post(i).Inter_p(2),~,stats] = ttest(data_no_smo_dbm,data_no_smo_vbm);
        post(i).Inter_t(2) = stats.tstat;
        
        % dbm: smo_no_smo
        [~,post(i).Inter_p(3),~,stats] = ttest(data_smo_dbm,data_no_smo_dbm);
        post(i).Inter_t(3) = stats.tstat;
        
        % vbm: smo_no_smo
        [~,post(i).Inter_p(4),~,stats] = ttest(data_smo_vbm,data_no_smo_vbm);
        post(i).Inter_t(4) = stats.tstat;
        
        post(i).Inter_fdr = gretna_FDR(post(i).Inter_p,p_thr);
        if ~isempty(post(i).Inter_fdr)
            % smo: dbm_vbm
            if post(i).Inter_p(1) <= post(i).Inter_fdr
                if post(i).Inter_t(1) > 0
                    post(i).Inter_pattern = [post(i).Inter_pattern 'Smo: DBM > VBM; '];
                else
                    post(i).Inter_pattern = [post(i).Inter_pattern 'Smo: DBM < VBM; '];
                end
            end
            
            % no_smo: dbm_vbm
            if post(i).Inter_p(2) <= post(i).Inter_fdr
                if post(i).Inter_t(2) > 0
                    post(i).Inter_pattern = [post(i).Inter_pattern 'No_smo: DBM > VBM; '];
                else
                    post(i).Inter_pattern = [post(i).Inter_pattern 'No_smo: DBM < VBM; '];
                end
            end
            
            % smo: dbm_vbm
            if post(i).Inter_p(3) <= post(i).Inter_fdr
                if post(i).Inter_t(3) > 0
                    post(i).Inter_pattern = [post(i).Inter_pattern 'DBM: Smo > No_smo; '];
                else
                    post(i).Inter_pattern = [post(i).Inter_pattern 'DBM: Smo < No_smo; '];
                end
            end
            
            % no_smo: dbm_vbm
            if post(i).Inter_p(4) <= post(i).Inter_fdr
                if post(i).Inter_t(4) > 0
                    post(i).Inter_pattern = [post(i).Inter_pattern 'VBM: Smo > No_smo; '];
                else
                    post(i).Inter_pattern = [post(i).Inter_pattern 'VBM: Smo < No_smo; '];
                end
            end
        end
    end
    
    data_smo_dbm = data.smo_dbm(net_mask);
    data_smo_vbm = data.smo_vbm(net_mask);
    
    data_no_smo_dbm = data.no_smo_dbm(net_mask);
    data_no_smo_vbm = data.no_smo_vbm(net_mask);
    
    post(i).data_smo_dbm = data_smo_dbm;
    post(i).data_smo_vbm = data_smo_vbm;
    post(i).data_no_smo_dbm = data_no_smo_dbm;
    post(i).data_no_smo_vbm = data_no_smo_vbm;
end

%% save
cd(path)
save('Connectivity_compare_ranova.mat','results','post')