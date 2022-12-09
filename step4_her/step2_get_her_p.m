%% get heritability and p-value of each edge
load('heritability_VBM.mat')
her_real = heritability(1,:);
her_rand = heritability(2:end,:);

roi_num = 48;
edge_mask = logical(tril(ones(roi_num),-1));

her = zeros(roi_num);
her(edge_mask) = her_real;
her = her + her';

her_p_vector = (sum(repmat(her_real,size(her_rand,1),1) <= her_rand)+1)/(size(her_rand,1)+1);
her_p = zeros(roi_num);
her_p(edge_mask) = her_p_vector;
her_p = her_p + her_p';

save('heritability_matrix_p_VBM.mat','her','her_p')