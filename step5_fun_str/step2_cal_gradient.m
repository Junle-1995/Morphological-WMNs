%% fun
load('unrelated_sub_matrix_fun.mat')
sub_mean_fun = mean(sub_matrix,3);

sub_cell_fun = cell(size(sub_matrix,3),1);
for i_sub = 1:length(sub_cell_fun)
    sub_cell_fun{i_sub} = sub_matrix(:,:,i_sub);
end

Gt_fun = GradientMaps('kernel','cosine','approach','dm');
Gt_fun = Gt_fun.fit(sub_mean_fun);

Gs_fun = GradientMaps('kernel','cosine','approach','dm','alignment','procrustes');
Gs_fun = Gs_fun.fit(sub_cell_fun','reference',Gt_fun.gradients{1});


%% save
save('unrelated_gradient.mat','Gt_fun','Gs_fun','Gt_str')