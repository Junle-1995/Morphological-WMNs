%% corr the connectivity and behav data via PLS
roi_num = 48;

path = pwd;

load('unrelated_sub_matrix_str_VBM.mat')
Mat_group = sub_matrix;
Mask_net = logical(tril(ones(roi_num),-1));
Mat_data = zeros(size(Mat_group,3),sum(Mask_net(:)));
for isub = 1:size(Mat_group,3)
    Mat_sub = Mat_group(:,:,isub);
    Mat_data(isub,:) = Mat_sub(Mask_net);
end
Mat_data = zscore(Mat_data);

rand_time = 10000;

BH_index = 1:6;

results = struct;

for i = BH_index
    cd(path)
    bev_data = load([pwd '\BH\BHtext' num2str(i) '.txt']);
    
    [~,~,XS,YS,~,PCTVAR,~,stats] = plsregress(Mat_data,zscore(bev_data),1);
    
    results(i).BH_index = i;
    results(i).XS = XS;
    results(i).YS = YS;
    results(i).variance_X = PCTVAR(1,:);
    results(i).variance_Y = PCTVAR(2,:);
    results(i).weighte = stats.W;
    
    for irand = 1:rand_time
        if mod(irand,floor(rand_time/10)) == 0
            disp(['Now calculating data in BH_' num2str(i) ' (' num2str(irand) '/' num2str(rand_time) ')  |' datestr(clock)])
        end
        
        rand_index = randperm(length(bev_data));
        while isequal(rand_index,1:length(bev_data)), rand_index = randperm(length(bev_data)); end
        
        bev_data_rand = bev_data(rand_index);
        [~,~,~,~,~,PCTVAR,~,~] = plsregress(Mat_data,zscore(bev_data_rand),1);
        results(i).variance_Y_rand(irand,1) = PCTVAR(2,:);
    end
    results(i).p = (sum(results(i).variance_Y_rand >= results(i).variance_Y)+1)/(rand_time+1);
end

cd(path)
save('Edge_BH_PLS_VBM.mat','results')