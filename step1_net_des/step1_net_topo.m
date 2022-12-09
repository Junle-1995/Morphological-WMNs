%% calculate the topological properties of morphological WM networks
path = pwd;
cd ..
data_path = [pwd '\step0_network_construction'];

cd(path)
[Net_para, Node_para] = gretna_batch_networkanalysis_bin([data_path '\VBM_smooth'], 'JS', 0.09, 0.3, 0.02, 100, 'pos', 's');
save('VBM_Net_para_bin_smooth.mat','Net_para','Node_para')

cd(path)
[Net_para, Node_para] = gretna_batch_networkanalysis_bin([data_path '\VBM_nosmooth'], 'JS', 0.09, 0.3, 0.02, 100, 'pos', 's');
save('VBM_Net_para_bin_nosmooth.mat','Net_para','Node_para')

cd(path)
[Net_para, Node_para] = gretna_batch_networkanalysis_bin([data_path '\DBM_smooth'], 'JS', 0.09, 0.3, 0.02, 100, 'pos', 's');
save('DBM_Net_para_bin_smooth.mat','Net_para','Node_para')

cd(path)
[Net_para, Node_para] = gretna_batch_networkanalysis_bin([data_path '\DBM_nosmooth'], 'JS', 0.09, 0.3, 0.02, 100, 'pos', 's');
save('DBM_Net_para_bin_nosmooth.mat','Net_para','Node_para')