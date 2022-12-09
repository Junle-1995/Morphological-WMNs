%% binary_classification to HC and NMO with edge
roi_num = 48;
load('group_data_DBM.mat')

edge_data = cat(3,data.HC,data.MS,data.NMO);

net_mask = logical(tril(ones(roi_num),-1));
Features = zeros(size(edge_data,3),sum(net_mask(:)));
for isub = 1:size(edge_data,3)
    edge_data_i = edge_data(:,:,isub);
    Features(isub,:) = edge_data_i(net_mask);
end

Group_label = [0*ones(size(data.HC,3),1);3*ones(size(data.MS,3),1);1*ones(size(data.NMO,3),1)];
Features = Features(Group_label ~= 3,:);
Group_label = Group_label(Group_label ~= 3);

K_fold = 10;
Method = 'SVM';
C_type = 'PCA';
C_thr = 80;

rep_time = 100;
rand_time = 1000;

results = struct;

for i_rep = 1:rep_time
    if mod(i_rep,floor(rep_time/10)) == 0
        disp(['Now calculating the data of repeation (' ...
            num2str(i_rep) '\' num2str(rep_time) ')  |' datestr(clock)])
    end
    
    Results = gretna_binary_classification_linear(Features, Group_label, K_fold, Method, C_type, C_thr);
    results.matrix(:,:,i_rep) = Results.Matrix;
    results.weight(:,i_rep) = Results.Consensus_weights;
    results.acc(:,i_rep) = Results.Accuracy;
    results.sen(:,i_rep) = Results.Sensitivity;
    results.spe(:,i_rep) = Results.Specificity;
    results.Explained_variance(i_rep,1) = mean(Results.Explained_variance);
end
results.Matrix = mean(results.matrix,3);
results.acc = mean(results.acc);
results.sen = mean(results.sen);
results.spe = mean(results.spe);


for i_rand = 1:rand_time
    if mod(i_rand,floor(rand_time/10)) == 0
        disp(['Now calculating the data of randmization (' ...
            num2str(i_rand) '\' num2str(rand_time) ')  |' datestr(clock)])
    end
    
    rand_index = randperm(length(Group_label));
    while isequal(rand_index,1:length(Group_label)), rand_index = randperm(length(Group_label)); end
    Group_label_rand = Group_label(rand_index);
    
    Results = gretna_binary_classification_linear(Features, Group_label_rand, K_fold, Method, C_type, C_thr);
    results.matrix_rand(:,:,i_rand) = Results.Matrix;
    results.acc_rand(i_rand,1) = Results.Accuracy;
    results.sen_rand(i_rand,1) = Results.Sensitivity;
    results.spe_rand(i_rand,1) = Results.Specificity;
end
results.acc_p = (sum(results.acc_rand >= results.acc)+1)/(rand_time+1);
results.sen_p = (sum(results.sen_rand >= results.sen)+1)/(rand_time+1);
results.spe_p = (sum(results.spe_rand >= results.spe)+1)/(rand_time+1);

save(['Edge_cla_HC_NMO_' Method '_' C_type '_' num2str(C_thr) '_DBM.mat'],'results')