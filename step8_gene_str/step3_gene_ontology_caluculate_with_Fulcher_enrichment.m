%% para
%% run after matlab 2018 version
file_name_all = {'Edge_GLM_gene_str_dis_commu_gene_PLS_DBM.mat',...
    'Edge_GLM_gene_str_dis_commu_gene_PLS_VBM.mat'};

for ifile = 1:length(file_name_all)
    file_name = file_name_all{ifile};
    load(file_name)

    z_thr = -icdf('normal',0.05,0,1);
    GO_size_min = 0;
    sizeFilter = [10 1000];
    load('gene_name_ID.mat')
    geneEntrezIDs = gene_name_ID(~isnan(gene_name_ID));


    %% get GO table
    GOTable = GetFilteredGOData('human-direct','biological_process',sizeFilter,geneEntrezIDs);


    %% get the gene index to each GO
    ID_index = struct;
    for iGO = 1:size(GOTable,1)
        ID_index(iGO).GOlabel = GOTable.GOIDlabel{iGO};
        ID_index(iGO).GOname = GOTable.GOName{iGO};
        ID_index(iGO).GOsize = GOTable.size(iGO);

        ID_index(iGO).index = find(ismember(geneEntrezIDs,GOTable.annotations{iGO}));
    end
    N = length(geneEntrezIDs);


    %% calculate the real GO score
    GO_results_pos = struct;
    GO_results_neg = struct;
    geneScores = zscore(gen_results.gen_weight);
    geneScores = geneScores(~isnan(gene_name_ID));
    gene_pos = find(geneScores >= z_thr);
    gene_neg = find(geneScores <= -z_thr);

    disp(['Now calculating the real GO score in   |' datestr(clock)])
    for iGO = 1:size(GOTable,1)
        GO_results_pos(iGO).GOlabel = GOTable.GOIDlabel{iGO};
        GO_results_pos(iGO).GOname = GOTable.GOName{iGO};

        GO_results_neg(iGO).GOlabel = GOTable.GOIDlabel{iGO};
        GO_results_neg(iGO).GOname = GOTable.GOName{iGO};

        GO_results_pos(iGO).GOsize = sum(ismember(ID_index(iGO).index,gene_pos));
        GO_results_neg(iGO).GOsize = sum(ismember(ID_index(iGO).index,gene_neg));
        GO_results_pos(iGO).GOscore = (N/length(ID_index(iGO).index))/(length(gene_pos)/GO_results_pos(iGO).GOsize);
        GO_results_neg(iGO).GOscore = (N/length(ID_index(iGO).index))/(length(gene_neg)/GO_results_neg(iGO).GOsize);

        if isnan(GO_results_pos(iGO).GOscore), GO_results_pos(iGO).GOscore = 0; end
        if isnan(GO_results_neg(iGO).GOscore), GO_results_neg(iGO).GOscore = 0; end
    end


    %% get the permutated GO score
    per_time = size(gen_results.gen_weight_moran,2);
    for iper = 1:per_time
        if mod(iper,floor(per_time/100)) == 0
            disp(['Now calculating the random GO score (' num2str(iper) '/' num2str(per_time) ')   |' datestr(clock)])
        end

        if corr(gen_results.gen_weight,gen_results.gen_weight_moran(:,iper)) < 0
            gen_results.gen_weight_moran(:,iper) = -gen_results.gen_weight_moran(:,iper);
        end
        geneScores = zscore(gen_results.gen_weight_moran(:,iper));
        geneScores = geneScores(~isnan(gene_name_ID));

        gene_pos = find(geneScores >= z_thr);
        gene_neg = find(geneScores <= -z_thr);

        for iGO = 1:size(GOTable,1)
            GOscore_pos_rand = (N/length(ID_index(iGO).index))/(length(gene_pos)/sum(ismember(ID_index(iGO).index,gene_pos)));
            GOscore_neg_rand = (N/length(ID_index(iGO).index))/(length(gene_neg)/sum(ismember(ID_index(iGO).index,gene_neg)));

            if isnan(GOscore_pos_rand), GOscore_pos_rand = 0; end
            if isnan(GOscore_neg_rand), GOscore_neg_rand = 0; end

            GO_results_pos(iGO).GOscore_rand(iper,1) = GOscore_pos_rand;
            GO_results_neg(iGO).GOscore_rand(iper,1) = GOscore_neg_rand;
        end
    end


    %% calculated the p-value of each GO
    for iGO = 1:size(GOTable,1)
        GO_results_pos(iGO).GOscore_per_p = mean(GO_results_pos(iGO).GOscore_rand >= GO_results_pos(iGO).GOscore);
        GO_results_neg(iGO).GOscore_per_p = mean(GO_results_neg(iGO).GOscore_rand >= GO_results_neg(iGO).GOscore);

        GO_results_pos(iGO).GOscore_Z_p = 1 - ...
            normcdf(GO_results_pos(iGO).GOscore,mean(GO_results_pos(iGO).GOscore_rand),std(GO_results_pos(iGO).GOscore_rand));
        GO_results_neg(iGO).GOscore_Z_p = 1 - ...
            normcdf(GO_results_neg(iGO).GOscore,mean(GO_results_neg(iGO).GOscore_rand),std(GO_results_neg(iGO).GOscore_rand));
    end

    % pos
    all_size = cat(2,GO_results_pos(:).GOsize);
    GO_results_pos(all_size < GO_size_min) = [];

    per_p_fdr = mafdr(cat(2,GO_results_pos(:).GOscore_per_p),'BHFDR',true,'showPlot',false);
    z_p_fdr = mafdr(cat(2,GO_results_pos(:).GOscore_Z_p),'BHFDR',true,'showPlot',false);
    [~,change_index_per] = sort(cat(2,GO_results_pos(:).GOscore_per_p),'ascend');

    for iGO = 1:length(GO_results_pos)
        GO_results_pos(iGO).GOscore_per_p_FDR = per_p_fdr(iGO);
        GO_results_pos(iGO).GOscore_Z_p_FDR = z_p_fdr(iGO);
    end
    GO_results_pos = GO_results_pos(change_index_per);

    % neg
    all_size = cat(2,GO_results_neg(:).GOsize);
    GO_results_neg(all_size < GO_size_min) = [];

    per_p_fdr = mafdr(cat(2,GO_results_neg(:).GOscore_per_p),'BHFDR',true,'showPlot',false);
    z_p_fdr = mafdr(cat(2,GO_results_neg(:).GOscore_Z_p),'BHFDR',true,'showPlot',false);
    [~,change_index_per] = sort(cat(2,GO_results_neg(:).GOscore_per_p),'ascend');

    for iGO = 1:length(GO_results_neg)
        GO_results_neg(iGO).GOscore_per_p_FDR = per_p_fdr(iGO);
        GO_results_neg(iGO).GOscore_Z_p_FDR = z_p_fdr(iGO);
    end
    GO_results_neg = GO_results_neg(change_index_per);


    %% save
    save([file_name(1:end-4) '_gene_ontology_pos_neg_enrichment.mat'],'GO_results_pos','GO_results_neg')
end