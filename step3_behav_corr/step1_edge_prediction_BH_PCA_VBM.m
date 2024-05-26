%% corr the connectivity and behav data via gretna_prediction
roi_num = 48;
load('unrelated_sub_matrix_str_VBM.mat')
Mask_net = logical(tril(ones(roi_num),-1));
Predictors = zeros(size(sub_matrix,3),sum(Mask_net(:)));
for isub = 1:size(sub_matrix,3)
    matrix_i = sub_matrix(:,:,isub);
    Predictors(isub,:) = matrix_i(Mask_net);
end

K_fold = 10;
C_type = 'PCA';
C_thr = 80;
R_type = 'glm';
Ncomp = [];
Lambda = [];

rep_time = 100;
rand_time = 10000;

path = pwd;
BH_path = [pwd '\BH'];
BH_index = 1:6;

results = struct;

for iBH = BH_index
    results(iBH).BH_index = iBH;
    
    cd(path)
    cd ..
    Responses = load([BH_path '\BHtext' num2str(iBH) '.txt']);
    
    for i_rep = 1:rep_time
        if mod(i_rep,floor(rep_time/10)) == 0
            disp(['Now calculating the data of repeation (' ...
                num2str(i_rep) '\' num2str(rep_time) ')  in ' num2str(iBH) ' |' datestr(clock)])
        end
        
        Results = gretna_prediction(Predictors, Responses, K_fold, C_type, C_thr, R_type, Ncomp, Lambda);
        results(iBH).Explained_variance(i_rep,1) = mean(Results.Explained_variance);
        results(iBH).R(i_rep,1) = Results.R;
        results(iBH).P(i_rep,1) = Results.P;
        results(iBH).Consensus_weights(:,i_rep) = Results.Consensus_weights;
    end
    
    for irand = 1:rand_time
        if mod(irand,floor(rand_time/10)) == 0
            disp(['Now calculating the data of randomization (' ...
                num2str(irand) '\' num2str(rand_time) ')  in ' num2str(iBH) ' |' datestr(clock)])
        end
        
        rand_index =randperm(length(Responses));
        while isequal(rand_index,1:length(Responses)), rand_index =randperm(length(Responses)); end
        
        Responses_rand = Responses(rand_index);
        
        Results = gretna_prediction(Predictors, Responses_rand, K_fold, C_type, C_thr, R_type, Ncomp, Lambda);
        results(iBH).R_rand(irand,1) = Results.R;
    end
    
    results(iBH).P_rand = (sum(results(iBH).R_rand >= mean(results(iBH).R)) + 1)/(rand_time + 1);
    
    cd(path)
    save(['Edge_BH_prediction_PCA_' num2str(C_thr) '_VBM.mat'],'results')
end