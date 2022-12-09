%% degree distribution
load('Sub_topo_str_VBM.mat')
roi_num = 48;

Glo_topo = {'acp','alp','acpratio','alpratio','aQ','aQratio'};
itopo = 1;
net_type = {'bin_no_smo','bin_smo'};

Nbins = 10;

results = struct;
save_path = [pwd '\pic'];
mkdir(save_path)

% degree
topo_data = [mean(sub_topo_bin_no_smo)',mean(sub_topo_bin_smo)'];
topo_data = topo_data(length(Glo_topo)+(itopo-1)*roi_num+1:length(Glo_topo)+itopo*roi_num,:);

for i = 1:length(net_type)
    results(i).net = net_type{i};
    [results(i).Para, results(i).R2] = gretna_degree_distribution(topo_data(:,i), Nbins);

    print([save_path '\Degree_dis_' net_type{i} '_VBM.tiff'],'-dtiff','-r600',gcf)
    close(gcf)
end

save('Degree_dis_VBM.mat','results')